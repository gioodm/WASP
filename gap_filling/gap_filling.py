#!/usr/bin/env python3.10
# @author Giorgia Del Missier


import argparse
import re, ast
import numpy as np
import networkx as nx
import itertools

np.random.seed(0)

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True,
                    help='foldseek results (.m8 file)')
parser.add_argument('--input_db', required=True,
                    help='foldseek results (.m8 file) containing tmscores for allvsall db search')
parser.add_argument('--input_rxn2up', required=True,
                    help='file containing rxn to uniprot mapping of orphan reactions')
parser.add_argument('--output', required=True,
                    help='output file name (.txt) containing rxn best hits in target taxid')
parser.add_argument('--eval_thr', required=True,
                    help='evalue threshold')
parser.add_argument('--bits_thr', required=True,
                    help='bitscore threshold')
parser.add_argument('--tm_thr', required=True,
                    help='tmscore threshold')
args = parser.parse_args()


def parse_fs(infile, eval_thr, bits_thr, tm_thr):
    """
    parse_fs() parses the foldseek output files
    """

    all_queries = dict()

    with open(infile) as fs_out:
        for line in fs_out:
            line = line.strip().split()
            query, target, evalue, bitscore, tmscore = line[0], line[1], float(line[-2]), int(line[-1]), float(line[-3])

            # extract query and target IDs from headers
            query = (lambda query : query.split("-")[1] if "-" in query else query)(query)
            target = (lambda target : target.split("-")[1] if "-" in target else target)(target)

            if query != target:
                if bitscore > bits_thr and evalue < eval_thr and tmscore > tm_thr:
                    try:
                        all_queries[query].append((target, tmscore))
                    except KeyError:
                        all_queries[query] = [(target, tmscore)]

    return all_queries


def fill_gaps(allq, rxn2up, tms, fout):
    '''
    fill_gaps() takes as input the structural alignment results produced by foldseek between the target organism proteome and
    known proteins associated with specific reactions (identified by Rhea IDs/EC numbers). Using this information, the function
    identifies hits in the target organisms which are structurally similar to the previously associated ones and can be used 
    to "fill-in" the gaps in GEM models.
    '''

    # define the column headers for the output file
    cols = ["#rxn ID", "rxn (extended name)", "rxn codes (Rhea ID/EC number)", "UniProt IDs (other organisms)", "top hit (target organism)", "other hits (target organism)"]
    
    # open the output file for writing and write the column headers
    with open(fout, 'w') as f_out:
        f_out.write("\t".join(cols) + "\n")

        # loop over all the Rhea IDs/EC numbers in rhea2up dictionary
        for key, values in rxn2up.items():

            # if there are no UniProt IDs associated with this Rhea ID/EC number or no Rhea ID/EC number
            # associated with this UniProt ID, write empty fields to output file
            if values[0] == [] or values[1] == []:
                fields = [key[0], key[1], '', '', '', '']
                f_out.write("\t".join(str(item) for item in fields) + "\n")

            # if there is only one UniProt ID associated with this Rhea ID/EC number
            elif len(values[1]) == 1:
                rxn_codes = ', '.join([str(x) for x in values[0] if x])
                # if the UniProt ID has hits in the target organism, sort them by similarity score and write to output file
                if values[1][0] in allq:
                    hits = sorted(allq[values[1][0]], key=lambda x: x[1], reverse=True)
                    fields = [key[0], key[1], rxn_codes, values[1][0], hits[0], str(hits[1:]).strip('[]')]
                    f_out.write("\t".join(str(item) for item in fields) + "\n")
                # if the UniProt ID does not have hits in the target organism, write empty fields to output file
                else:
                    fields = [key[0], key[1], rxn_codes, values[1][0], '', '']
                    f_out.write("\t".join(str(item) for item in fields) + "\n")

            else:
                rxn_codes = ', '.join([str(x) for x in values[0] if x])

                # create a graph containing subclusters of structurally similar proteins within all the 
                # UniProt IDs matching a certain reaction code
                G = nx.Graph()
                comb = list(itertools.combinations(values[1], 2))
                for e in comb:
                    if e in tms:
                        G.add_edge(e[0], e[1])

                # get the connected components in the graph
                Gcc = sorted(nx.connected_components(G), key=len, reverse=True)

                # if there are no connected components, write out the result to the output file
                if Gcc == []:
                    up = ', '.join([str(x) for x in values[1] if x])
                    fields = [key[0], key[1], rxn_codes, up, '', '']
                    f_out.write("\t".join(str(item) for item in fields) + "\n")
                else:
                    for e in Gcc:
                        # get the query hits for each node in the component
                        hits = [allq[j] for j in e if j in allq]
                        # if there are no hits, write out the result to the output file
                        if hits == []:
                            up = ', '.join([str(x) for x in values[1] if x])
                            fields = [key[0], key[1], rxn_codes, up, '', '']
                            f_out.write("\t".join(str(item) for item in fields) + "\n")
                        # otherwise, calculate the average similarity score for each ID and write out the result to the output file
                        else:
                            id_dict = dict()
                            for h in hits:
                                for tup in h:
                                    if tup[0] not in id_dict:
                                        id_dict[tup[0]] = [tup[1]]
                                    else:
                                        id_dict[tup[0]].append(tup[1])
                            output = sorted([(k, round(sum(v)/len(v),3)) for k, v in id_dict.items() if len(v) == len(hits)], key=lambda x: x[1], reverse=True)
                            up = ', '.join([str(x) for x in e if x])
                            if output != []:
                                fields = [key[0], key[1], rxn_codes, up, output[0], str(output[1:]).strip('[]')]
                                f_out.write("\t".join(str(item) for item in fields) + "\n")
                            else:
                                fields = [key[0], key[1], rxn_codes, up, '', '']
                                f_out.write("\t".join(str(item) for item in fields) + "\n")


# parsing the input files
all_queries = parse_fs(args.input, float(args.eval_thr), int(args.bits_thr), float(args.tm_thr))

rxn2up_dict = dict()
with open(args.input_rxn2up) as f_rxn2up:
    for line in f_rxn2up:
        if line.startswith("#"):
            pass
        else:
            line = line.strip().split("\t")
            rxn_codes = [item for item in ast.literal_eval(line[2]) if item not in ['', '[]']]
            rxn_up = [item for item in ast.literal_eval(line[3]) if item not in ['', '[]']]
            rxn2up_dict[(line[0], line[1])] = (rxn_codes, rxn_up)

with open(args.input_db) as f_db:
    tms_dict = {(line.strip().split("\t")[0].split("-")[1], line.strip().split("\t")[1].split("-")[1]):float(line.strip().split("\t")[-3]) \
            for line in f_db if float(line.strip().split("\t")[-3]) > float(args.tm_thr)}

fill_gaps(all_queries, rxn2up_dict, tms_dict, args.output)



