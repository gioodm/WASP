#!/usr/bin/env python3.11
# @author Giorgia Del Missier

#%%

import requests, re
from requests.adapters import HTTPAdapter, Retry

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
    url = f"{UNIPROT_API}/uniprotkb/search?query=accession={pnames}&fields=accession,gene_oln,protein_name&format=tsv&size=500"

    for batch, total in get_batch(url):
        for r in batch.text.splitlines()[1:]:
            r = r.split('\t')
            annotation[r[0]] = r[1:]

    return annotation

#%%

import pandas as pd
import numpy as np

gem_out = pd.read_csv("results/1111708/1111708_hits.txt", sep='\t', header=0)

gem_filtered = gem_out.dropna(subset=["rxn codes (Rhea ID/EC number)"])
gem_filtered = gem_out.dropna(subset=["top hit (target organism)"])

tot = 0
fout = open("results/1111708_hits_new.txt", "w")
for index, row in gem_filtered.iterrows():
    proteins, tmscores = list(), dict()
    top = eval(row["top hit (target organism)"])
    proteins.append(top[0])
    tmscores[top[0]] = (top[1])
    if type(row["other hits (target organism)"]) != float and len(row["other hits (target organism)"].split("), ")) > 1:
        other = [i[0] for i in eval(row["other hits (target organism)"])]
        other_tm = [i[1] for i in eval(row["other hits (target organism)"])]
        proteins.extend(other)
        for i, j in enumerate(other):
            tmscores[j] = other_tm[i]
    elif type(row["other hits (target organism)"]) != float and len(row["other hits (target organism)"].split("), ")) == 1:
        other = eval(row["other hits (target organism)"])
        proteins.append(other[0])
        tmscores[other[0]] = other[1]

    tot += len(proteins)
    results = get_UniProt(proteins)

    for r in proteins:
        new_row = [row['#rxn ID'], row['rxn (extended name)'], row['rxn codes (Rhea ID/EC number)'], row['UniProt IDs (other organisms)'], r, tmscores[r], results[r][0], results[r][1]]
        fout.write("\t".join(str(item) for item in new_row) + "\n")

fout.close()

# %%
