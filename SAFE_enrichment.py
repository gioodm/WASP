#!/usr/bin/env python3.10
#@author Giorgia Del Missier

import argparse, sys, os
import pandas as pd
import numpy as np
import multiprocessing
import altair as alt
from altair_saver import save

np.random.seed(0)

sys.path.append(f'{os.getcwd()}/safepy')
import safe

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True,
                    help='input: file containing complete annotation for each cluster (.txt file)')
parser.add_argument('--radius', required=True,
                    help='neighborhood radius for SAFE analysis')
parser.add_argument('--edgelist', required=True,
                    help='input: file containing network edge list with respective weights (.txt file)')
parser.add_argument('--wd', required=True,
                    help='working directory')
parser.add_argument('--output', required=True,
                    help='output: file name containing improved annotation')
parser.add_argument('--outfig', required=True,
                    help='output: figure name for summary statistics of improved annotation')
parser.add_argument('--nan', required=False,
                    help='output: file name containing IDs of uncharacterised proteins')
parser.add_argument('--identifiers', required=True, 
                    help="foldseek file to retrieve UniProt identifiers from target organism")
parser.add_argument('--iteration', required=True, 
                    help="iteration number")
args = parser.parse_args()


def process_identifier(args_list):
	i, df = args_list
	subdf = df[['#Cluster', 'UniProt ID', i]]
	flattened_ids = subdf[i].dropna().astype(str).str.split(';').explode()
	final_matrix = pd.DataFrame()
	if not flattened_ids.empty:
		unique_ids = flattened_ids.dropna().unique()
		final_matrix = pd.DataFrame(columns=unique_ids)

		for index, row in subdf.iterrows():
			protein = row['UniProt ID']
			ids = row[i]
			if pd.isna(ids):
				ids = list()
			final_matrix.loc[protein] = [1 if j in ids else 0 for j in unique_ids]

	output_file = f"{args.wd}/SAFE/{i}_matrix.txt"
	final_matrix.to_csv(output_file, sep='\t')


def perform_safe_analysis(i, wd, radius, sf, enrichment_thr = 0.05):
	sf.load_attributes(attribute_file=f'SAFE/{i}_matrix.txt')
	sf.define_neighborhoods(neighborhood_radius=float(radius))
	sf.compute_pvalues(background='network', num_permutations=1000, processes=64)
	sf.print_output_files(output_dir=f'{wd}/SAFE/{i}_')
    
	df_tmp = pd.read_csv(f"{wd}/SAFE/{i}_node_properties_annotation.txt", sep="\t")
        
	significant = dict()
	for index, row in df_tmp.iterrows():
		key = row['key']
		significant[key] = [(col, round(row[col],3)) for col in df_tmp.columns[3:] if row[col] > -np.log10(enrichment_thr) and not col.startswith("Unnamed")]
    
	return i, significant


def prepare_df(d):

	df = pd.DataFrame()

	for i in d:
		tmp_df = pd.DataFrame(d[i]).transpose()
		tmp_df = tmp_df.stack().reset_index()
		tmp_df.columns = ['id', 'group', 'count']
		tmp_df['category'] = i

		df = pd.concat([df, tmp_df])

	return df


def make_barcharts(txid, coll):

	df_taxid = prepare_df(txid)
	df_collateral = prepare_df(coll)

	chart_taxid = alt.Chart(df_taxid).mark_bar().encode(

		# tell Altair which field to group columns on
    	x=alt.X('group:N', title=None, sort=None),  # Add sort argument

    	# tell Altair which field to use as Y values and how to calculate
    	y=alt.Y('sum(count):Q',
        	axis=alt.Axis(
            	grid=False,
            	title=None)),

    	# tell Altair which field to use to use as the set of columns to be represented in each group
    	column=alt.Column('id:N', sort=alt.SortField(field='id', order='ascending'), title=None),  # Add sort argument

    	# tell Altair which field to use for color segmentation 
    	color=alt.Color('category:N', sort=None,
            scale=alt.Scale(
                # make it look pretty with an enjoyable color pallet
                range=['#001219','#005F73','#94D2BD','#EE9B00','#BB3E03','#9B2226'],
            )),

    	order=alt.Order(
      	# Sort the segments of the bars by this field
      	'category',
      	sort='ascending')
        	)\
        	.configure_view(
        	# remove grid lines around column clusters
        	strokeOpacity=0    
    	).properties(
    	width=150)

	chart_taxid = chart_taxid.configure_legend(labelFontSize=8)
	chart_taxid.save(f"{args.outfig}_taxid_iter{args.iteration}.png", ppi=200)

	chart_collateral = alt.Chart(df_collateral).mark_bar().encode(

		# tell Altair which field to group columns on
    	x=alt.X('group:N', title=None, sort=None),  # Add sort argument

    	# tell Altair which field to use as Y values and how to calculate
    	y=alt.Y('sum(count):Q',
        	axis=alt.Axis(
            	grid=False,
            	title=None)),

    	# tell Altair which field to use to use as the set of columns to be represented in each group
    	column=alt.Column('id:N', sort=alt.SortField(field='id', order='ascending'), title=None),  # Add sort argument

    	# tell Altair which field to use for color segmentation 
    	color=alt.Color('category:N', sort=None,
            scale=alt.Scale(
                # make it look pretty with an enjoyable color pallet
                range=['#001219','#005F73','#94D2BD','#EE9B00','#BB3E03','#9B2226'],
            )),

    	order=alt.Order(
      	# Sort the segments of the bars by this field
      	'category',
      	sort='ascending')
        	)\
        	.configure_view(
        	# remove grid lines around column clusters
        	strokeOpacity=0    
    	).properties(
    	width=150)

	chart_collateral = chart_collateral.configure_legend(labelFontSize=8)
	chart_collateral.save(f"{args.outfig}_collateral_iter{args.iteration}.png", ppi=200)


def analyze_safe_results(df_anno, df_safe, ids, allp):

	taxid, collateral = {}, {}

	for key in ["total", "previously annotated", "not annotated", "nan2nan", "new annotation by WASP", "nan2nan in dark clusters"]:
		taxid[key], collateral[key] = {}, {}

		for identifier in ids:
			taxid[key][f'{ids.index(identifier) + 1}.{identifier}'] = {'total': 0, 'annotation': 0, 'nan2nan / new': 0, 'nan2nan (dark)': 0}
			collateral[key][f'{ids.index(identifier) + 1}.{identifier}'] = {'total': 0, 'annotation': 0, 'nan2nan / new': 0, 'nan2nan (dark)': 0}

	list_nan2nan = list()

	cols = ["UniProt ID"] + ids + ["Organism"]
	taxid_outdf, collateral_outdf = pd.DataFrame(columns=cols), pd.DataFrame(columns=cols)

	nclusters = df_anno["#Cluster"].unique()
	for c in nclusters:
		
		subdf_anno = df_anno[df_anno["#Cluster"] == c]
		subdf_safe = df_safe[df_safe["UniProt ID"].isin(subdf_anno["UniProt ID"])]

		dark_anno = subdf_anno[ids].apply(lambda col: col.isna().all()).to_dict()
		dark_safe = subdf_safe[ids].apply(lambda col: col.isna().all()).to_dict()
		
		tmp_taxid, tmp_collateral = pd.DataFrame(columns=cols), pd.DataFrame(columns=cols)

		for i in ids:

			taxid['total'][f'{ids.index(i) + 1}.{i}']['total'] += len(subdf_anno[subdf_anno['UniProt ID'].isin(allp['UniProt ID'])])
			collateral['total'][f'{ids.index(i) + 1}.{i}']['total'] += len(subdf_anno[~subdf_anno['UniProt ID'].isin(allp['UniProt ID'])])

			if dark_anno[i] and dark_safe[i]:
				proteins = subdf_anno[["UniProt ID", "Organism"]]
				proteins_taxid = proteins[proteins['UniProt ID'].isin(allp['UniProt ID'])]
				list_nan2nan.extend(proteins_taxid["UniProt ID"].tolist())

				taxid['nan2nan'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += len(proteins_taxid)
				taxid['not annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += len(proteins_taxid)
				taxid['nan2nan in dark clusters'][f'{ids.index(i) + 1}.{i}']['nan2nan (dark)'] += len(proteins_taxid)

				collateral['nan2nan'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += len(proteins) - len(proteins_taxid)
				collateral['not annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += len(proteins) - len(proteins_taxid)
				collateral['nan2nan in dark clusters'][f'{ids.index(i) + 1}.{i}']['nan2nan (dark)'] += len(proteins) - len(proteins_taxid)

			else:

				pre = subdf_anno[['UniProt ID', i, "Organism"]].reset_index(drop=True)
				new = subdf_safe[['UniProt ID', i]].reset_index(drop=True)
				new = new.sort_values(by='UniProt ID', key=lambda x: x.map(dict(zip(pre['UniProt ID'], pre.index)))).reset_index(drop=True)

				pre_taxid = pre[pre['UniProt ID'].isin(allp['UniProt ID'])]

				# compute new annotation by WASP
				changes = pre[i].isna() & new[i].notna()
				changed_rows = new[changes]
				# compute nan2nan
				changes_n2n = pre[i].isna() & new[i].isna()
				unchanged_rows_n2n = new[changes_n2n]

				changed_rows = changed_rows.merge(pre[changes][['UniProt ID', 'Organism']], on='UniProt ID', how='left')
				changed_rows_taxid = changed_rows[changed_rows['UniProt ID'].isin(allp['UniProt ID'])]
				changed_rows_collateral = changed_rows[~changed_rows['UniProt ID'].isin(allp['UniProt ID'])]

				unchanged_rows_n2n = unchanged_rows_n2n.merge(pre[changes_n2n][['UniProt ID', 'Organism']], on='UniProt ID', how='left')
				unchanged_rows_taxid_n2n = unchanged_rows_n2n[unchanged_rows_n2n['UniProt ID'].isin(allp['UniProt ID'])]
				unchanged_rows_collateral_n2n = unchanged_rows_n2n[~unchanged_rows_n2n['UniProt ID'].isin(allp['UniProt ID'])]
				list_nan2nan.extend(unchanged_rows_taxid_n2n["UniProt ID"].tolist())

				taxid['not annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += (len(changed_rows_taxid) + len(unchanged_rows_taxid_n2n))
				taxid['new annotation by WASP'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += len(changed_rows_taxid)
				taxid['nan2nan'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += len(unchanged_rows_taxid_n2n)
				taxid['previously annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += len(pre_taxid) - (len(changed_rows_taxid) + len(unchanged_rows_taxid_n2n))

				collateral['not annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += (len(changed_rows_collateral) + len(unchanged_rows_collateral_n2n))
				collateral['new annotation by WASP'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += len(changed_rows_collateral)
				collateral['nan2nan'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += len(unchanged_rows_collateral_n2n)
				collateral['previously annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += (len(pre) - len(pre_taxid)) - (len(changed_rows_collateral) + len(unchanged_rows_collateral_n2n))

				tmp_taxid = pd.concat([tmp_taxid, changed_rows_taxid], ignore_index=True, sort=False)
				tmp_collateral = pd.concat([tmp_collateral, changed_rows_collateral], ignore_index=True, sort=False)

		tmp_taxid = tmp_taxid.fillna('')
		tmp_collateral = tmp_collateral.fillna('')

		tmp_taxid = tmp_taxid.groupby('UniProt ID').agg({
                    'Pfam': ''.join,
                    'PANTHER': ''.join,
                    'CATH': ''.join,
                    'EC number': ''.join,
                    'Rhea ID': ''.join,
                    'GO terms': ''.join,
                    'Organism': 'first'
                    }).reset_index()

		tmp_collateral = tmp_collateral.groupby('UniProt ID').agg({
                    'Pfam': ''.join,
                    'PANTHER': ''.join,
                    'CATH': ''.join,
                    'EC number': ''.join,
                    'Rhea ID': ''.join,
                    'GO terms': ''.join,
                    'Organism': 'first'
                    }).reset_index()

		taxid_outdf = pd.concat([taxid_outdf, tmp_taxid], axis=0, ignore_index=True)
		collateral_outdf = pd.concat([collateral_outdf, tmp_collateral], axis=0, ignore_index=True)

	make_barcharts(taxid, collateral)

	return taxid_outdf, collateral_outdf, list_nan2nan


if __name__ == '__main__':

	identifiers = ["Pfam", "PANTHER", "CATH", "EC number", "Rhea ID", "GO terms"] 
	anno_df = pd.read_csv(args.input, sep="\t", header=0)

	args_list = [(i, anno_df) for i in identifiers]

    # Use multiprocessing.Pool to run the loop in parallel
	with multiprocessing.Pool() as pool:
		pool.map(process_identifier, args_list)

	# Initialize the safe object
	sf = safe.SAFE(path_to_safe_data=f'{args.wd}/')

	# Load and save the network
	sf.load_network(network_file=f'{args.edgelist}')

	with multiprocessing.Pool() as pool:
		results = pool.starmap(perform_safe_analysis, [(i, args.wd, args.radius, sf) for i in identifiers])
        
	safe_enriched = dict(results)
	safe_df = pd.DataFrame.from_dict(safe_enriched, orient='index').transpose()
	for i in identifiers:
		safe_df[i] = safe_df[i].apply(lambda x: sorted(x, key=lambda x: x[1], reverse=True))
		safe_df[i] = safe_df[i].apply(lambda x: ';'.join(f'({item[0]}, {item[1]})' for item in x) if len(x) > 0 else np.nan)
	safe_df = safe_df.reset_index().rename(columns={'index': 'UniProt ID'})

	allproteins = list()
	with open(args.identifiers) as fup:
		for line in fup:
			line = line.strip().split()
			query = line[0]
			query = (lambda query : query.split("-")[1] if "-" in query else query)(query)
			query = (lambda query : query.split(".gz")[0][:-4] if ".gz" in query else query)(query)
			allproteins.append(query)

	allproteins = list(set(allproteins))  # Ensure unique entries
	allproteins = pd.DataFrame({'UniProt ID': allproteins})

	outdf_taxid, outdf_collateral, nan2nan = analyze_safe_results(anno_df, safe_df, identifiers, allproteins)

	try:
		with pd.ExcelWriter(f'{args.output}_taxid.xlsx', engine='openpyxl', mode='a') as writer:
			outdf_taxid.to_excel(writer, sheet_name=f'Iteration_{args.iteration}', index=False)
		with pd.ExcelWriter(f'{args.output}_collateral.xlsx', engine='openpyxl', mode='a') as writer:
			outdf_collateral.to_excel(writer, sheet_name=f'Iteration_{args.iteration}', index=False)
	except:
		with pd.ExcelWriter(f'{args.output}_taxid.xlsx', engine='xlsxwriter') as writer:
			outdf_taxid.to_excel(writer, sheet_name=f'Iteration_{args.iteration}', index=False)
		with pd.ExcelWriter(f'{args.output}_collateral.xlsx', engine='xlsxwriter') as writer:
			outdf_collateral.to_excel(writer, sheet_name=f'Iteration_{args.iteration}', index=False)

	print(f"{len(set(nan2nan))} could not be annotated using SAFE... trying with increased number of neighbors")
	
	# write nan2nan IDs to a .txt file
	with open(args.nan2nan_out, "a") as fnan:
		for i in set(nan2nan):
			fnan.write(i + "\n")
	fnan.close()
