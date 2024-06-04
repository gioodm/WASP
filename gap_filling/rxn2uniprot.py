#!/usr/bin/env python3.11
# @author Giorgia Del Missier


import argparse
import requests


parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True,
                    help='.txt file containing orphan rxn IDs, rxn extended names and rxn code (Rhea ID/EC number)')
parser.add_argument('--output', required=True,
                    help='output file name (.txt) containing mapping of orphan rxns to corresponding UniProt IDs')
parser.add_argument('--output_ids', required=True,
                    help='output file name (.txt) containing AlphaFold IDs of all UniProt IDs mapped')
args = parser.parse_args()


headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_5) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/50.0.2661.102 Safari/537.36'}


def get_UniProt(protein):
    """
    get_UniProt() takes a RheaID/EC number identifier as input and uses the UniProt REST API to retrieve all the UniProt IDs that are associated
    with that catalytic function/reaction ID.
    """

    UNIPROT_API = "https://rest.uniprot.org/uniprotkb/"

    try:
        # If the input protein ID starts with "rhea", it is a Rhea ID
        if protein.startswith("RHEA"):
            rhea_id = protein.split(":")[1]
            # query UniProt API for UniProt IDs associated with the catalytic activity of the given Rhea ID
            uniprot_response = requests.get(f"{UNIPROT_API}stream?format=list&query=%28%28cc_catalytic_activity_exp%3A%22rhea%3A{rhea_id}%22%29%20AND%20%28database%3Aalphafolddb%29%29", headers=headers)

        # if the input protein ID starts with "EC", it is an EC number
        elif protein.startswith("EC"):
            ec_id = protein.split(":")[1]
            # query UniProt API for UniProt IDs associated with the catalytic activity of the given EC number
            uniprot_response = requests.get(f"{UNIPROT_API}stream?format=list&query=%28%28cc_catalytic_activity_exp%3A%22EC%3A{ec_id}%22%29%20AND%20%28database%3Aalphafolddb%29%29", headers=headers)

        # extract the UniProt IDs from the response
        uniprot_ids = uniprot_response.text.strip().split()
        return uniprot_ids

    except:
        # Return an empty string if there is an error
        return ''


# define list to store all UniProt IDs retrieved
all_uniprot_ids = []

with open(args.output_ids, "w") as f_ids:
    with open(args.output, "w") as f_out:
        cols = ["#rxn ID", "rxn (extended name)", "rxn codes (Rhea ID/EC number)", "UniProt IDs"]
        f_out.write("\t".join(cols) + "\n")
        rhea2uniprot = {}

        with open(args.input) as f_orphans2rhea:
            for line in f_orphans2rhea:
                if line.startswith("#"):
                    # skip comment lines
                    pass
                else:
                    # split line into orphan and associated rxn codes
                    line = line.strip().split("\t")
                    if len(line) > 1:
                        orphan, rxn_name, rxn_codes = line[0], line[1], line[2:]
                    else:
                        orphan, rxn_name, rxn_codes = line[0], '', []

                    # initialize list to store UniProt IDs associated with this orphan
                    orphan_uniprot_ids = []

                    # iterate over rxn codes associated with this orphan
                    for rxn in rxn_codes:
                        rxn = rxn.strip(" , ").split(", ")
                        for r in rxn:
                            # if the UniProt IDs associated with this rxn code were previously retrieved, use them from the dictionary
                            if r in rhea2uniprot:
                                uniprot_ids = rhea2uniprot[r]
                            # otherwise, retrieve
                            else:
                                uniprot_ids = get_UniProt(r)
                                # store UniProt IDs in the dictionary for future use
                                rhea2uniprot[r] = uniprot_ids
                                # add the retrieved UniProt IDs to the list of UniProt IDs associated with this orphan
                                orphan_uniprot_ids.extend(uniprot_ids)

                    # remove duplicates
                    orphan_uniprot_ids = list(set(orphan_uniprot_ids))
                    # add to the list of all UniProt IDs
                    all_uniprot_ids.extend(orphan_uniprot_ids)

                    # write to the output file
                    fields = [orphan, rxn_name] + [rxn_codes] + [orphan_uniprot_ids]
                    f_out.write("\t".join(str(item) for item in fields) + "\n")

            # remove duplicates from the list of all UniProt IDs retrieved
            all_uniprot_ids = list(set(all_uniprot_ids))
            # write to the output file
            for uniprot_id in all_uniprot_ids:
                f_ids.write(f"AF-{uniprot_id}-F1-model_v4\n")

