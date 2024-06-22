#!/usr/bin/env python3.10
# @Author: Giorgia Del Missier

import argparse
import requests


parser = argparse.ArgumentParser(description="Map orphan reactions to UniProt IDs using Rhea or EC numbers.")
parser.add_argument('--input', required=True, help='.txt file containing orphan rxn IDs, rxn extended names, and rxn code (Rhea ID/EC number)')
parser.add_argument('--output', required=True, help='Output file name (.txt) containing mapping of orphan rxns to corresponding UniProt IDs')
parser.add_argument('--output_ids', required=True, help='Output file name (.txt) containing AlphaFold IDs of all UniProt IDs mapped')
args = parser.parse_args()


def get_UniProt(protein):
    """
    Retrieve UniProt IDs associated with a given Rhea ID or EC number using the UniProt REST API.
    """
    UNIPROT_API = "https://rest.uniprot.org/uniprotkb/"

    try:

        if protein.startswith("RHEA"):
            rhea_id = protein.split(":")[1]
            # Query UniProt API for UniProt IDs associated with the Rhea ID
            response = requests.get(f"{UNIPROT_API}stream?format=list&query=(cc_catalytic_activity_exp:\"rhea:{rhea_id}\") AND (database:alphafolddb)")

        elif protein.startswith("EC"):
            ec_id = protein.split(":")[1]
            # Query UniProt API for UniProt IDs associated with the EC number
            response = requests.get(f"{UNIPROT_API}stream?format=list&query=(cc_catalytic_activity_exp:\"EC:{ec_id}\") AND (database:alphafolddb)")

        # Extract UniProt IDs from the response
        uniprot_ids = response.text.strip().split()
        return uniprot_ids
    
    except Exception as e:
        print(f"Error retrieving UniProt IDs for {protein}: {e}")
        return []


all_uniprot_ids = []


with open(args.output_ids, "w") as f_ids, open(args.output, "w") as f_out:

    cols = ["#rxn ID", "rxn (extended name)", "rxn codes (Rhea ID/EC number)", "UniProt IDs"]
    f_out.write("\t".join(cols) + "\n")
    
    # Dictionary to cache UniProt IDs for Rhea/EC numbers
    rhea2uniprot = {}

    # Read the input file line by line
    with open(args.input) as f_orphans2rhea:
        for line in f_orphans2rhea:
            if line.startswith("#"):
                # Skip comment lines
                continue

            # Split line into components
            line = line.strip().split("\t")
            orphan = line[0]
            rxn_name = line[1] if len(line) > 1 else ''
            rxn_codes = line[2:] if len(line) > 2 else []

            # List to store UniProt IDs for the current orphan reaction
            orphan_uniprot_ids = []

            # Process each reaction code associated with the orphan reaction
            for rxn in rxn_codes:
                rxn = rxn.strip(" , ").split(", ")
                for r in rxn:
                    # Use cached UniProt IDs if available
                    if r in rhea2uniprot:
                        uniprot_ids = rhea2uniprot[r]
                    else:
                        # Retrieve UniProt IDs and cache them
                        uniprot_ids = get_UniProt(r)
                        rhea2uniprot[r] = uniprot_ids
                    
                    # Add retrieved UniProt IDs to the orphan reaction's list
                    orphan_uniprot_ids.extend(uniprot_ids)

            # Remove duplicates
            orphan_uniprot_ids = list(set(orphan_uniprot_ids))
            # Add to the global list of all UniProt IDs
            all_uniprot_ids.extend(orphan_uniprot_ids)

            # Write the orphan reaction and its associated UniProt IDs to the output file
            fields = [orphan, rxn_name, ", ".join(rxn_codes), ", ".join(orphan_uniprot_ids)]
            f_out.write("\t".join(fields) + "\n")

    # Remove duplicates from the global list of UniProt IDs
    all_uniprot_ids = list(set(all_uniprot_ids))

    # Write AlphaFold IDs to the output IDs file
    for uniprot_id in all_uniprot_ids:
        f_ids.write(f"AF-{uniprot_id}-F1-model_v4\n")
