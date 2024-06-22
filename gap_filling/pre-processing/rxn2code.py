#!/usr/bin/env python3.11
# @Author: Giorgia Del Missier

import argparse
import pandas as pd
import requests


# Constants for API URLs
RHEA_API = "https://www.rhea-db.org/rhea/"
BIGG_API = "http://bigg.ucsd.edu/api/v2/universal/reactions/"

def parse_arguments():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(description='Map orphan reactions to Rhea and EC numbers.')
    parser.add_argument('--input', required=True, help='.tsv input file with orphan reactions and annotations (divided into cols)')
    parser.add_argument('--cols', required=True, help='Column IDs in the annotation file')
    parser.add_argument('--mxn_xref', required=True, help='.tsv input file with MetaNetX code external references')
    parser.add_argument('--output', required=True, help='Output file name (.txt) containing orphan reactions mapped to Rhea and EC numbers')
    return parser.parse_args()

def request_rhea_data(query):
    """
    Request Rhea data based on a query string.
    """
    try:
        response = requests.get(f"{RHEA_API}?query={query}&columns=rhea-id&format=list")
        return response.text.strip().split("\n")
    except requests.RequestException:
        return []

def request_ec_data(query):
    """
    Request EC number data based on a query string.
    """
    try:
        response = requests.get(f"{RHEA_API}?query={query}&columns=enzyme-class&format=tsv")
        ec_data = response.text.strip().split("\n")[1:]  # Skip the header row
        return [e.split(" ")[0] for e in ec_data]  # Extract only the EC numbers
    except requests.RequestException:
        return []

def process_kegg(keggid):
    """
    Process KEGG ID to retrieve corresponding Rhea IDs and EC numbers.
    """
    rhea = request_rhea_data(f"KEGG%3A{keggid}")
    ec = request_ec_data(f"KEGG%3A{keggid}")
    return rhea, ec

def process_bigg(biggid):
    """
    Process BiGG ID to retrieve corresponding Rhea IDs and EC numbers.
    """
    all_rhea, all_ec = [], []
    try:
        response = requests.get(f"{BIGG_API}{biggid}")
        payload = response.json()
        
        # Check if there are database links in the response
        if "database_links" in payload:
            if "RHEA" in payload["database_links"]:
                # Retrieve Rhea IDs and corresponding EC numbers
                for link in payload["database_links"]["RHEA"]:
                    r = link["id"]
                    all_rhea.append(f"RHEA:{r}")
                    ec = request_ec_data(f"RHEA%3A{r}")
                    all_ec.extend(ec)

            if "KEGG Reaction" in payload["database_links"]:
                # Retrieve Rhea IDs and EC numbers from KEGG links
                for link in payload["database_links"]["KEGG Reaction"]:
                    keggid = link["id"]
                    rhea, ec = process_kegg(keggid)
                    all_rhea.extend(rhea)
                    all_ec.extend(ec)
    except requests.RequestException:
        pass
    return list(set(all_rhea)), list(set(all_ec))

def process_pubmed(pmedid):
    """
    Process PubMed ID to retrieve corresponding Rhea IDs and EC numbers.
    """
    rhea = request_rhea_data(f"pubmed%3A{pmedid}")
    ec = request_ec_data(f"pubmed%3A{pmedid}")
    return rhea, ec

def parse_mxn_xref(file):
    """
    Parse MetaNetX cross-reference file to create a dictionary.
    """
    mxn_dict = {}
    with open(file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            # Check for relevant prefixes and store in dictionary
            if fields[0].startswith(('biggR:', 'rheaR:', 'keggR:')):
                key = fields[1]
                value = fields[0]
                if key in mxn_dict:
                    mxn_dict[key].append(value)
                else:
                    mxn_dict[key] = [value]
    return mxn_dict

def process_metanetx(mxnid, mxn_dict):
    """
    Process MetaNetX ID to retrieve corresponding Rhea IDs and EC numbers.
    """
    all_rhea, all_ec = [], []
    if mxnid in mxn_dict:
        for entry in mxn_dict[mxnid]:
            rhea, ec = [], []
            if entry.startswith('biggR:'):
                rhea, ec = process_bigg(entry.split(":")[1])
            elif entry.startswith('keggR:'):
                rhea, ec = process_kegg(entry.split(":")[1])
            elif entry.startswith("rheaR:"):
                rhea = [f"RHEA:{entry.split(':')[1]}"]
                ec = request_ec_data(f"RHEA%3A{entry.split(':')[1]}")
            all_rhea.extend(rhea)
            all_ec.extend(ec)
    return list(set(all_rhea)), list(set(all_ec))

def process_columns(filename, columns, mxn_dict):
    """
    Process the input file to map reactions to Rhea IDs and EC numbers.
    """
    cols = ["#rxn", "rxn (extended name)", "Rhea ID", "EC number"]
    df = pd.read_csv(filename, sep='\t')
    newdf = pd.DataFrame(columns=cols)

    for _, row in df.iterrows():
        all_rheas, all_ecs = [], []
        for col in columns:
            val = row[col]
            if isinstance(val, str) and val.strip():
                # Handle cases where the cell contains a list-like string
                vals = val[1:-1].split(", ") if val.startswith('[') else [val]
                vals = [v.replace("'", "") for v in vals]
                for v in vals:
                    rhea, ec = [], []
                    if col == 'kegg.reaction':
                        rhea, ec = process_kegg(v)
                    elif col == 'bigg.reaction':
                        rhea, ec = process_bigg(v)
                    elif col == 'pubmed':
                        rhea, ec = process_pubmed(v)
                    elif col == 'metanetx.reaction':
                        rhea, ec = process_metanetx(v, mxn_dict)
                    elif col == "ec-code":
                        ec = [f"EC:{v}"]
                    elif col == "rhea":
                        rhea = [f"RHEA:{v}"]

                    all_rheas.extend(rhea)
                    all_ecs.extend(ec)

        # Create a new row with the results
        new_line = {
            "#rxn": row["rxn_id"],
            "rxn (extended name)": row["rxn_name"],
            "Rhea ID": ', '.join(set(all_rheas)),
            "EC number": ', '.join(set(all_ecs))
        }
        newdf = pd.concat([newdf, pd.DataFrame([new_line], columns=cols)], ignore_index=True)
    return newdf

def main():
    """
    Main function to execute the script.
    """
    args = parse_arguments()
    columns = args.cols.split(", ")
    mxn_dict = parse_mxn_xref(args.mxn_xref)
    df = process_columns(args.input, columns, mxn_dict)
    df.to_csv(args.output, sep='\t', index=False)
    print(f"Processed data written to {args.output}")

if __name__ == '__main__':
    main()
