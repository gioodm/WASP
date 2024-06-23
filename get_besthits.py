#!/usr/bin/env python3.10
# Author: Giorgia Del Missier

import argparse

# Argument parser setup
parser = argparse.ArgumentParser(description="Extract best hits from Foldseek results based on evalue and bitscore thresholds.")
parser.add_argument('--input', required=True, help='Input Foldseek results file (.m8 format)')
parser.add_argument('--output', required=True, help='Output file to store IDs of best hits')
parser.add_argument('--eval_thr', required=True, type=float, help='Evalue threshold for filtering hits')
parser.add_argument('--bits_thr', required=True, type=int, help='Bitscore threshold for filtering hits')
parser.add_argument('--hit_pos', required=False, type=int, default=1, help='Hit position to extract (1st, 2nd, etc.)')
args = parser.parse_args()

def get_besthit(input_file, n, evalue_threshold, bitscore_threshold):
    """
    Extract the best hit for each query based on evalue and bitscore thresholds.
    
    Parameters:
    - input_file: Path to the input .m8 file
    - n: The nth hit to extract for each query
    - evalue_threshold: The maximum evalue for a hit to be considered
    - bitscore_threshold: The minimum bitscore for a hit to be considered
    
    Returns:
    - A set of the best hit targets
    """
    queries = dict()

    with open(input_file) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            query, target, evalue, bitscore = line[0], line[1], float(line[-2]), int(line[-1])
            
            # Strip potential suffix from query and target identifiers
            query_stripped = query.split("-")[1] if "-" in query else query
            target_stripped = target.split("-")[1] if "-" in target else target
            
            # Ensure query and target are not the same and apply thresholds
            if query_stripped != target_stripped:
                if evalue <= evalue_threshold and bitscore >= bitscore_threshold: 
                    # Store targets for each query
                    if query not in queries:
                        queries[query] = []
                    if target not in queries[query]:
                        queries[query].append(target)
    
    # Extract the nth hit for each query if it exists
    hits = [targets[n-1] for query, targets in queries.items() if len(targets) >= n]
    print(f"\nHits found for {len(hits)} queries\n")

    return set(hits)

# Get the best hits using the provided arguments
best_hits = get_besthit(args.input, args.hit_pos, args.eval_thr, args.bits_thr)

# Write the best hits to the output file
with open(args.output, "w") as f_out:
    for hit in best_hits:
        f_out.write(hit + "\n")

print(f"Best hits have been written to {args.output}")
