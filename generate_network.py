#!/usr/bin/env python3.10
# @author Giorgia Del Missier


import argparse
import json
import numpy as np
import networkx as nx

# set random seed for reproducibility
np.random.seed(0)


parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True,
                    help='input: foldseek parsed results (.json file)')
parser.add_argument('--input_bh', required=True,
                    help='input: foldseek parsed reciprocal search results (.json file)')
parser.add_argument('--nan', required=True,
                    help='input: IDs of uncharacterised proteinsfrom previous step (.txt file)')
parser.add_argument('--neighbours', required=True,
                    help='maximum number of neighbours included in the network for each query')
parser.add_argument('--output', required=True,
                    help='output: .txt file containing clusters and additional metrics')
parser.add_argument('--edgelist', required=True,
                    help='output: .txt file containing network edge list')
args = parser.parse_args()


def create_network(all_queries, reciprocal_queries, max_hits):
    """
    create_network() creates a network graph containing clusters of proteins with similar structure
    """

    G = nx.Graph()

    # add edges for each query and its top n neighbours
    for query, hits in all_queries.items():

        if hits[0][0] in reciprocal_queries:
            reci_hits = {e[0] : e[2] for e in reciprocal_queries[hits[0][0]][:max_hits]}
            if query in reci_hits:
                G.add_edge(query, hits[0][0], weight=hits[0][2])

                i = 1
                while i < len(hits) and i < max_hits:
                    if hits[i][0] in reci_hits:
                        G.add_edge(query, hits[i][0], weight=hits[i][2])
                        G.add_edge(hits[0][0], hits[i][0], weight=reci_hits[hits[i][0]])
                    i+=1

    return G


if __name__ == "__main__":

    with open(args.input) as fall, open(args.input_bh) as freci:
        all_queries = json.load(fall)
        reciprocal_queries = json.load(freci)

    with open(args.nan) as fnan:
        nan2nan = [line.strip() for line in fnan]

    if nan2nan != []:
        all_queries = {key: value for key, value in all_queries.items() if key in nan2nan}
        
    G = create_network(all_queries, reciprocal_queries, int(args.neighbours))

    diff = set(all_queries.keys()).difference(set(G.nodes()))

    # write nan2nan IDs to a .txt file
    with open(args.nan, "w") as foutnan:
        for i in diff:
            foutnan.write(i + "\n")
    foutnan.close()

    print(f"Found RBSH hits for {len(all_queries) - len(diff)} out of {len(all_queries)} proteins from the previous step\n")
    print(f"{len(diff)} proteins in the target organism had no RBSH hits... trying again with increased number of neighbours\n")

    print(f"Network statistics generated using {len(all_queries) - len(diff)} RBSH hits:")
    print(f"Number of nodes: {len(G.nodes())}")
    print(f"Number of edges: {len(G.edges())}\n")

    print(f"Number of generated clusters: {nx.number_connected_components(G)}")

    # sort the clusters by size (number of nodes)
    clusters_sorted = [c for c in sorted(nx.connected_components(G), key=len, reverse=True)]

    # write the clusters information to output file
    columns = ["#Cluster", "Number of Nodes", "Number of Edges", "Average Degree", "Average Clustering Coefficient", "Node List"]
    output_file = open(args.output, "w")
    output_edgelist = open(args.edgelist, "w")
    output_file.write("\t".join(columns) + "\n")

    for idx, cluster in enumerate(clusters_sorted):

        S = G.subgraph(cluster)
        avg_degree = round(sum(dict(S.degree(weight='weight')).values()) / len(S), 3)
        avg_clustering_coef = round(nx.average_clustering(S, weight='weight'), 3)

        for edge in S.edges.data("weight"):
            output_edgelist.write("\t".join(str(item) for item in edge) + "\n")

        fields = [idx, len(S.nodes()), len(S.edges()), avg_degree, avg_clustering_coef, str(cluster).strip("{ }")]
        output_file.write("\t".join(str(item) for item in fields) + "\n")

    output_file.close()
    output_edgelist.close()


