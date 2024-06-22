#!/usr/bin/env python3.10
# @author Giorgia Del Missier

import argparse
import requests, re
import multiprocessing
import time


parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True,
                    help='.txt file containing communities information')
parser.add_argument('--output', required=True,
                    help='output file name (.txt) containing complete annotation for each hit')
args = parser.parse_args()


def get_UniProt(protein):
    """
    get_UniProt() takes a UniProt identifier as input and uses the UniProt REST API to retrieve information about the protein's function and annotation.
    It then processes the response and returns it as a list.
    """
    
    try:
        UNIPROT_API = "https://rest.uniprot.org/"

        # send a query to UniProt for the given protein
        uniprot = requests.get(f"{UNIPROT_API}/uniprotkb/search?query={protein}&fields=id,rhea,cc_function,annotation_score,organism_name&format=tsv", headers=headers)

        # extract the requested fields from the query response and split into a list
        fields = uniprot.text.strip('\n').split('\t')[5:]

        # if the protein has a "FUNCTION" annotation, extract it from the field
        if len(fields) == 4 and fields[1] != '':
            fields[1] = fields[1].split("FUNCTION: ")[1]

    except:
        fields = ['', '', '', '']

    # return the list of extracted and processed fields
    return fields


def process_line(line):

    result_list = []
    
    if line.startswith("#"):
        # skip comment lines
        return None
    else:
        fields = line.strip().split("\t")

        # if there is only one hits
        if len(fields) == 5:
            rxn = fields[0]
            rxn_name = fields[1]
            rxn_code = fields[2]
            top_hit = fields[4].strip("()").split(",")[0].strip("''")

            # get annotation
            uniprot_fields = get_UniProt(top_hit)

            # combine all the fields and add them to the result list
            fields = [rxn, rxn_name, rxn_code, top_hit, pfam, panther, ko] + uniprot_fields
            result_list.append("\t".join(str(field) for field in fields) + "\n")
            time.sleep(1)

        # if there are multiple hits
        elif len(fields) == 6:
            rxn = fields[0]
            rxn_name = fields[1]
            rxn_code = fields[2]
            top_hit = fields[4].strip("()").split(",")[0].strip("''")
            other_hits = fields[5]
            pattern = r"\('(\w+)',"
            other_hits = re.findall(pattern, other_hits)
            hits = [top_hit] + other_hits
            for h in hits:
                # get annotation
                pfam, panther = get_InterPro(h)
                ko = get_KO(h)
                uniprot_fields = get_UniProt(h)
            
                # combine all the fields and add them to the result list
                fields = [rxn, rxn_name, rxn_code, h, pfam, panther, ko] + uniprot_fields
                result_list.append("\t".join(str(field) for field in fields) + "\n")
                time.sleep(1)
    
    return result_list


if __name__ == "__main__":
    
    # define the column headers for the output file
    cols = ["rxn ID", "rxn (extended name)", "rxn code (Rhea ID/EC number)", "UniProt ID", "Pfam", "Panther", "Kegg Orthology", "Rhea ID", "Function", "Annotation score", "Organism"]

    with open(args.input) as fin, open(args.output, "w") as fout:
        fout.write("\t".join(cols) + "\n")
        # read the lines of the input file
        lines = fin.readlines()
        
        # use multiprocessing to process the lines in parallel
        with multiprocessing.Pool(processes=16) as pool:
            results = pool.map(process_line, lines)
            for result in results:
                if result is not None and len(results) > 0:
                    # write the results to the output file
                    fout.write("".join(result))

