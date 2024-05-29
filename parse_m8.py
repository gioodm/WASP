#!/usr/bin/env python3.10
# @author Giorgia Del Missier


import argparse
import json


parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True,
                    help='input: foldseek results (.m8 file path)')
parser.add_argument('--input_bh', required=True,
                    help='input: foldseek reciprocal search results (.m8 file path)')
parser.add_argument('--input_normalisation', required=True,
                    help='input: foldseek results for bitscore normalisation (.m8 file path)')
parser.add_argument('--input_normalisation_bh', required=True,
                    help='input: foldseek best reciprocals results for bitscore normalisation (.m8 file path)')
parser.add_argument('--proteome_size', required=True, 
                    help="proteome size (i.e. number of AlphaFold models) of selected taxid")
parser.add_argument('--eval_thr', required=True,
                    help='evalue threshold')
parser.add_argument('--bits_thr', required=True,
                    help='bitscore threshold')
args = parser.parse_args()


def parse_fs(fin, fnormalisation, eval_thr, bits_thr):
    """
    parse_fs() parses the foldseek output files
    """

    queries, normalisation = dict(), dict()

    with open(fnormalisation) as fnorm:
        for line in fnorm:
            line = line.strip().split()
            query, target, evalue, bitscore = line[0], line[1], float(line[-2]), int(line[-1])
            
            # extract query and target IDs from headers
            query = (lambda query : query.split("-")[1] if "-" in query else query)(query)
            query = (lambda query : query.split(".gz")[0][:-4] if ".gz" in query else query)(query)
            target = (lambda target : target.split("-")[1] if "-" in target else target)(target)
            target = (lambda target : target.split(".gz")[0][:-4] if ".gz" in target else target)(target)

            if query == target:
                if query not in normalisation:
                    normalisation[query] = bitscore


    with open(fin) as fs_in:
        for line in fs_in:
            line = line.strip().split()
            query, target, evalue, bitscore = line[0], line[1], float(line[-2]), int(line[-1])
            
            # extract query and target IDs from headers
            query = (lambda query : query.split("-")[1] if "-" in query else query)(query)
            query = (lambda query : query.split(".gz")[0][:-4] if ".gz" in query else query)(query)
            target = (lambda target : target.split("-")[1] if "-" in target else target)(target)
            target = (lambda target : target.split(".gz")[0][:-4] if ".gz" in target else target)(target)
            
            if query != target:
                if bitscore > bits_thr and evalue < eval_thr:
                    bitscore = round(bitscore/normalisation[query], 3)
                    try:
                        if (target, evalue, bitscore) not in queries[query]:
                            queries[query].append((target, evalue, bitscore))
                    except KeyError:
                        queries[query] = [(target, evalue, bitscore)]

    return queries


if __name__ == "__main__":
    
    # parsing the input files
    all_queries = parse_fs(args.input, args.input_normalisation, float(args.eval_thr), int(args.bits_thr))
    reciprocal_queries = parse_fs(args.input_bh, args.input_normalisation_bh, float(args.eval_thr), int(args.bits_thr))

    print(f"Found significant hits for {len(all_queries)} out of {args.proteome_size} proteins in the target organism\n")

    with open(f"{args.input[:-3]}.json", "w") as fall, open(f"{args.input_bh[:-3]}.json", "w") as freci:
        json.dump(all_queries, fall)
        json.dump(reciprocal_queries, freci)