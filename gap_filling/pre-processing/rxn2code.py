#!/usr/bin/env python3.11
# @author Giorgia Del Missier


import argparse
import pandas as pd
import requests


parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True,
                    help='.txt input file with orphan reactions and annotations (divided into cols)')
parser.add_argument('--cols', required=True,
                    help='columns ids in the annotation file')
parser.add_argument('--mxn_xref', required=True,
                    help='.tsv input file with MetaNetX code external references \
                          retrieved from: https://www.metanetx.org/mnxdoc/mnxref.html, reac_xref.tsv')
parser.add_argument('--output', required=True,
                    help='output file name (.txt) containing orphan reactions mapped to rxn code (Rhea and EC number)')
args = parser.parse_args()


RHEA_API = "https://www.rhea-db.org/rhea/"


def process_kegg(keggid):

  try:
    rhea = requests.get(f"{RHEA_API}?query=KEGG%3A{keggid}&columns=rhea-id&format=list")
    rhea = rhea.text.strip().split("\n")

  except:
    rhea = []

  try:
    ec = requests.get(f"{RHEA_API}?query=KEGG%3A{keggid}&columns=enzyme-class&format=tsv")
    ec = ec.text.strip("\n").split("\n")[1:]
    ec = [e.split(" ")[0] for e in ec]

  except:
    ec = []

  return rhea, ec


def process_bigg(biggid):

  all_rhea, all_EC = [], []

  try:
    req = requests.get(f'http://bigg.ucsd.edu/api/v2/universal/reactions/{biggid}')
    payload = req.json()

    if payload["database_links"] != {}:

      if payload["database_links"]["RHEA"]:
        for j in payload["database_links"]["RHEA"]:
          r = j["id"]
          all_rhea.append("RHEA:" + r)
          try:
            ec = requests.get(f"{RHEA_API}?query=RHEA%3A{r}&columns=enzyme-class&format=tsv")
            ec = ec.text.strip("\n").split("\n")[1:]
            ec = list(set([e.split(" ")[0] for e in ec]))
            if ec != []:
              all_EC += ec
            
          except:
            pass

      if payload["database_links"]["KEGG Reaction"]:
          for j in payload["database_links"]["KEGG Reaction"]:
            keggid = j["id"]
            rhea, ec = process_kegg(keggid)
            if rhea != []:
              all_rhea += rhea
            if ec != []:
              all_EC += ec

  except:
    pass

  return all_rhea, all_EC


def process_metanetx(mxnid):
	pass

def process_pubmed(pmedid):

  try:
    rhea = requests.get(f"{RHEA_API}?query=pubmed%3A{pmedid}&columns=rhea-id&format=list")
    rhea = rhea.text.strip().split("\n")

  except:
    rhea = []

  try:
    ec = requests.get(f"{RHEA_API}?query=pubmed%3A{pmedid}&columns=enzyme-class&format=tsv")
    ec = ec.text.strip("\n").split("\n")[1:]
    ec = [e.split(" ")[0] for e in ec]

  except:
    ec = []
  
  return rhea, ec


def parse_mxn_xref(file):

  d = {}
  with open(file) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if fields[0].startswith(('biggR:', 'rheaR:', 'keggR:')):
          try:
            d[fields[1]].append(fields[0])
          except:
            d[fields[1]] = [fields[0]]
  return d


def process_metanetx(mxnid, d):

  all_rhea, all_ec = [], []
  if mxnid in d:
    for i in d[mxnid]:
      if i.startswith('biggR:'):
        i = i.split(":")[1]
        rhea, ec = process_bigg(i)
      elif i.startswith('keggR:'):
        i = i.split(":")[1]
        rhea, ec = process_kegg(i)
      elif i.startswith("rheaR:"):
        i = i.split(":")[1]
        rhea = ["RHEA:" + i]
        ec = requests.get(f"{RHEA_API}?query=RHEA%3A{i}&columns=enzyme-class&format=tsv")
        ec = ec.text.strip("\n").split("\n")[1:]
        ec = list(set([e.split(" ")[0] for e in ec]))

      if rhea != []:
        all_rhea += rhea
      if ec != []:
        all_ec += ec

  return list(set(all_rhea)), list(set(all_ec))


def process_columns(filename, columns, mxn_dict):
    cols = ["#rxn", "rxn (extended name)", "Rhea ID", "EC number"]
    newdf = pd.DataFrame(columns = cols)
    df = pd.read_csv(filename, sep='\t')
    for _, row in df.iterrows():
        all_rheas, all_ecs = [], []
        for col in columns:
            val = row[col]
            if isinstance(val, str) and val.strip() != '':
              if val.startswith('['):
                vals = val[1:-1].split(", ")
                vals = [e.replace("'", "") for e in vals]
              else:
                vals = [val]
              for val in vals:
                rhea, EC = [], []

                if col == 'kegg.reaction':
                    rhea, EC = process_kegg(val)

                elif col == 'bigg.reaction':       #or col == 'rxn_id':
                    rhea, EC = process_bigg(val)

                elif col == 'pubmed':
                    rhea, EC = process_pubmed(val) 

                elif col == 'metanetx.reaction':
                    rhea, EC = process_metanetx(val, mxn_dict)

                elif col == "ec-code":
                    EC = ["EC:"+val]

                elif col == "rhea":
                    rhea = ["RHEA:"+val]

                if rhea != []:
                    all_rheas += rhea
                if EC != []:
                    all_ecs += EC

        all_rheas = '{}'.format(', '.join(set(all_rheas)))
        all_ecs = '{}'.format(', '.join(set(all_ecs)))
        row_data = [row["rxn_id"], row["rxn_name"], all_rheas, all_ecs]
        new_line = pd.DataFrame([row_data], columns= cols)
        newdf = pd.concat([newdf, new_line])

    return newdf


cols = args.cols.split(", ")
mxn_d = parse_mxn_xref(args.mxn_xref)
df = process_columns(args.input, cols, mxn_d)
df.to_csv(args.output, sep='\t', index=False)
