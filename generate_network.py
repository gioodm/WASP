#!/usr/bin/env python3.10
# Author: Giorgia Del Missier

import argparse
import json
import numpy as np
import networkx as nx

# Set random seed for reproducibility
np.random.seed(0)

# Argument parser setup
parser = argparse.ArgumentParser(description="Generate a network of protein clusters based on Foldseek results.")
parser.add_argument('--input', required=True, help='Input Foldseek parsed results (.json file)')
parser.add_argument('--input_bh', required=True, help='Input Foldseek parsed reciprocal search results (.json file)')
parser.add_argument('--nan', required=True, help='Input file containing IDs of uncharacterized proteins from previous step (.txt file)')
parser.add_argument('--neighbours', required=True, type=int, help='Maximum number of neighbours included in the network for each query')
parser.add_argument('--output', required=True, help='Output file containing clusters and additional metrics (.txt file)')
parser.add_argument('--edgelist', required=True, help='Output file containing network edge list (.txt file)')
args = parser.parse_args()

def create_network(all_queries, reciprocal_queries, max_hits):
    """
    Create a network graph containing clusters of proteins with similar structure.

    Parameters:
    - all_queries: Dictionary of queries and their hits from the input file
    - reciprocal_queries: Dictionary of reciprocal hits from the input_bh file
    - max_hits: Maximum number of neighbours to include for each query

    Returns:
    - A networkx Graph object representing the protein clusters
    """
    G = nx.Graph()

    # Add edges for each query and its top n neighbours
    for query, hits in all_queries.items():
        if hits[0][0] in reciprocal_queries:
            reci_hits = {e[0]: e[2] for e in reciprocal_queries[hits[0][0]][:max_hits]}
            if query in reci_hits:
                G.add_edge(query, hits[0][0], weight=hits[0][2])

                i = 1
                while i < len(hits) and i < max_hits:
                    if hits[i][0] in reci_hits:
                        G.add_edge(query, hits[i][0], weight=hits[i][2])
                        G.add_edge(hits[0][0], hits[i][0], weight=reci_hits[hits[i][0]])
                    i += 1

    return G

if __name__ == "__main__":
    # Load data from input files
    with open(args.input) as fall, open(args.input_bh) as freci:
        all_queries = json.load(fall)
        reciprocal_queries = json.load(freci)

    # Load IDs of uncharacterized proteins
    with open(args.nan) as fnan:
        nan2nan = [line.strip() for line in fnan]

    # Filter queries to include only uncharacterized proteins
    if nan2nan:
        all_queries = {key: value for key, value in all_queries.items() if key in nan2nan}

    # Create the network graph
    G = create_network(all_queries, reciprocal_queries, args.neighbours)

    # Find the difference between all queries and the graph nodes
    diff = set(all_queries.keys()).difference(set(G.nodes()))

    # Write nan2nan IDs to the .txt file
    with open(args.nan, "w") as foutnan:
        for i in diff:
            foutnan.write(i + "\n")

    print(f"Found RBSH hits for {len(all_queries) - len(diff)} out of {len(all_queries)} proteins from the previous step")
    print(f"{len(diff)} proteins in the target organism had no RBSH hits... trying again with increased number of neighbours")

    # Print network statistics
    print(f"Network statistics generated using {len(all_queries) - len(diff)} RBSH hits:")
    print(f"Number of nodes: {len(G.nodes())}")
    print(f"Number of edges: {len(G.edges())}")
    print(f"Number of generated clusters: {nx.number_connected_components(G)}")

    # Sort the clusters by size (number of nodes)
    clusters_sorted = [c for c in sorted(nx.connected_components(G), key=len, reverse=True)]

    # Write the clusters information to output files
    columns = ["#Cluster", "Number of Nodes", "Number of Edges", "Average Degree", "Average Clustering Coefficient", "Node List"]
    with open(args.output, "w") as output_file, open(args.edgelist, "w") as output_edgelist:
        output_file.write("\t".join(columns) + "\n")

        for idx, cluster in enumerate(clusters_sorted):
            S = G.subgraph(cluster)
            avg_degree = round(sum(dict(S.degree(weight='weight')).values()) / len(S), 3)
            avg_clustering_coef = round(nx.average_clustering(S, weight='weight'), 3)

            for edge in S.edges.data("weight"):
                output_edgelist.write("\t".join(str(item) for item in edge) + "\n")

            fields = [idx, len(S.nodes()), len(S.edges()), avg_degree, avg_clustering_coef, str(cluster).strip("{ }")]
            output_file.write("\t".join(str(item) for item in fields) + "\n")

    print(f"Cluster information has been written to {args.output}")
    print(f"Edge list has been written to {args.edgelist}")
