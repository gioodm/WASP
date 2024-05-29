#!/usr/bin/env python3.10
# @author Giorgia Del Missier

import argparse
import requests
import re
from requests.adapters import HTTPAdapter, Retry

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True,
                    help='input: file containing cluster information (.txt file)')
parser.add_argument('--output', required=True,
                    help='output: .txt file containing complete annotation for each member of the clusters')
args = parser.parse_args()


re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=25, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url, timeout=25)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)


def get_UniProt(chunk):

    UNIPROT_API = "https://rest.uniprot.org/"
    pnames = "%20OR%20".join(chunk)

    annotation = dict()
    url = f"{UNIPROT_API}/uniprotkb/search?query=accession={pnames}&fields=accession,xref_pfam,xref_panther,xref_gene3d,ec,rhea,go,length,annotation_score,organism_name&format=tsv&size=500"
    for batch, total in get_batch(url):
        for r in batch.text.splitlines()[1:]:
            r = r.split('\t')

            r = [r.strip(";") for r in r]
            r[4] = r[4].replace("; ", ";")
            r[5] = r[5].replace(" ", ";")
            r[6] = ';'.join(['GO:' + g for g in re.findall(r'GO:(\d+)', r[6])])
            if ':' in r[2]:
                ids = r[2].split(";")
                r[2] = ';'.join(set([pthr for pthr in ids if ':' in pthr]))

            annotation[r[0]] = r[1:]

    return annotation


if __name__ == "__main__":

    cols = ["#Cluster", "UniProt ID", "Pfam", "PANTHER", "CATH", "EC number", "Rhea ID", "GO terms", "Length", "Annotation score", "Organism"]

    with open(args.input) as fin, open(args.output, "w") as fout:
        fout.write("\t".join(cols) + "\n")
        lines = fin.readlines()
        
        for line in lines[1:]:
            line = line.strip().split("\t")
            n = line[0]
            nodes = line[-1].split(",")
            nodes = [node.strip(" ' ' ") for node in nodes]

            up_nodes = [n for n in nodes if re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', n) and len(n) in (6, 10)]
            chunksize = 150
            chunks = [up_nodes[x:x+chunksize] for x in range(0, len(up_nodes), chunksize)]

            results = dict()
            for chunk in chunks:
                results.update(get_UniProt(chunk))

            for node in nodes:
                if node in results:
                    fields = [n, node] + results[node]
                    fout.write("\t".join(str(item) for item in fields) + "\n")
                else:
                    fields = [n, node, "", "", "", "", "", "", "", "", ""]
                    fout.write("\t".join(str(item) for item in fields) + "\n")
                    