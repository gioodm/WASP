#!/usr/bin/env python3.10
# @author Giorgia Del Missier


import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True,
                    help='input: foldseek results (.m8 file path)')
parser.add_argument('--output', required=True,
                    help='output: .txt file containing IDs of best hits')
parser.add_argument('--eval_thr', required=True,
                    help='evalue threshold')
parser.add_argument('--bits_thr', required=True,
                    help='bitscore threshold')
parser.add_argument('--hit_pos', required=False, default=1,
                    help='hit to extract (1st, 2nd, 3rd etc)')
args = parser.parse_args()


def get_besthit(input_file, n, evalue_threshold, bitscore_threshold):
    """ 
    get_besthit() extracts the best hit for each query
    """

    queries = dict()

    with open(input_file) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            query, target, evalue, bitscore = line[0], line[1], float(line[-2]), int(line[-1])
            query_stripped = (lambda query : query.split("-")[1] if "-" in query else query)(query)
            target_stripped = (lambda target : target.split("-")[1] if "-" in target else target)(target)
            if query_stripped != target_stripped:
                if evalue <= evalue_threshold and bitscore >= bitscore_threshold: 
                    try:
                        if target not in queries[query]:
                            queries[query].append(target)
                    except KeyError:
                        queries[query] = [target]


    hits = [targets[n-1] for query, targets in queries.items() if len(targets) >= n]
    print("\nHits found for " + str(len(hits)) + " queries\n")

    return set(hits)


best_hits = get_besthit(args.input, int(args.hit_pos), float(args.eval_thr), int(args.bits_thr))

# write best hits to a .txt file
with open(args.output, "w") as f_out:
    for hit in best_hits:
        f_out.write(hit + "\n")

f_out.close()

