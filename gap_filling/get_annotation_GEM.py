#!/usr/bin/env python3.11
# @author Giorgia Del Missier


import argparse
import requests, re
import multiprocessing
import requests_cache
import time


parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True,
                    help='.txt file containing communities information')
parser.add_argument('--output', required=True,
                    help='output file name (.txt) containing complete annotation for each hit')
args = parser.parse_args()


headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36'}

requests_cache.install_cache('my_cache')


def get_InterPro(protein):
    """
    get_InterPro() takes a UniProt identifier as input and uses InterPro API to retrieve information about 
    the protein's domains (Pfam) and subfamilies (Panther). The function returns a list of Pfam domains 
    and a string containing the Panther subfamily information. If the requests fail, the function returns 
    empty lists for the Pfam domains and the Panther subfamily. 
    """

    # initialize variables
    pfam, panther = [], ""

    try:
        # construct the URL for the InterPro API request and send the request
        url = f"https://www.ebi.ac.uk:443/interpro/wwwapi//entry/all/protein/unreviewed/{protein}"
        req = requests.get(url, headers=headers)

        # if the API request is successful, parse the JSON payload
        payload = req.json()
        # loop through all the InterPro entries in the payload
        for item in payload["results"]:
            # check if the entry has any metadata
            if item["metadata"]:
                try:
                    # extract the Pfam domains from the metadata and add them to the pfam list
                    dom = item["metadata"]["member_databases"]["pfam"]
                    pfam = pfam + list(dom.keys())
                except:
                    # if there are no Pfam domains, pass and continue
                    pass
            # check if the entry has any protein information
            if item["proteins"]:
                try:
                    # extract the subfamily information from the protein and store it as a comma-separated string in panther
                    panther = item["proteins"][0]['entry_protein_locations'][0]["subfamily"]
                    panther = ", ".join(list(panther.values()))
                except:
                    # if there is no subfamily information, pass and continue
                    pass

    except:
        try:
            # if the initial API request fails, try again with reviewed proteins
            url = f"https://www.ebi.ac.uk:443/interpro/wwwapi//entry/all/protein/reviewed/{protein}"
            req = requests.get(url, headers=headers)
            payload = req.json()

            for item in payload["results"]:
                if item["metadata"]:
                    try:
                        dom = item["metadata"]["member_databases"]["pfam"]
                        pfam = pfam + list(dom.keys())
                    except:
                        pass
                if item["proteins"]:
                    try:
                        panther = item["proteins"][0]['entry_protein_locations'][0]["subfamily"]
                        panther = ", ".join(list(panther.values()))
                    except:
                        pass
        except:
            # if both requests fail, return empty lists and strings
            pass

    return pfam, panther


def get_KO(protein):
    """
    get_KO() takes a UniProt identifier as input and retrieves the KEGG Orthology (KO) ID using the KEGG API.
    If no KO ID is successfully retrieved, an empty string is returned.
    """

    try:
        # build URL for KEGG API request with protein ID
        url = f"https://rest.kegg.jp/conv/genes/uniprot:{protein}"
    
        # send request to KEGG API
        req = requests.get(url, headers=headers)
        payload = req.text.strip().split()
    
        # if protein ID is valid, use it to retrieve KO (KEGG Orthology) ID
        if payload != [] and req.status_code == 200:
            # build new URL for KEGG API request with retrieved gene ID
            url = f"https://rest.kegg.jp/link/ko/{payload[1]}"
            req = requests.get(url)
            payload = req.text.strip().split()
            if payload != [] and req.status_code == 200:
                ko = payload[1].split(":")[1]
                return ko
    except:
        # return empty string if KO ID could not be retrieved
        return ""

    return ""


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
            pfam, panther = get_InterPro(top_hit)
            ko = get_KO(top_hit)
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

