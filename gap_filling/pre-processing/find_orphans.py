#!/usr/bin/env python3.11
# @author Giorgia Del Missier


import argparse
import cobra
import numpy as np
import pandas as pd
import ast


parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True,
                    help='GEM model')
parser.add_argument('--extension', required=True,
                    help='GEM model extension (xml, sbml, mat, json)')
parser.add_argument('--output', required=True,
                    help='output file name (.txt) containing orphan reactions+annotation')
args = parser.parse_args()


def find_orphan_rxns(model):
    '''
    find_orphan_rxns() identifies all orphan reaction in the GEM model (i.e. reactions with no associated genes),
    excluding exchange reactions.
    '''

    # create an empty list to store the IDs of orphan reactions
    orphan_reactions = []
    index = 1

    # iterate over all reactions in the model
    for reaction in model.reactions:
        # check if the reaction has no associated genes
        if len(reaction.genes) == 0:
            # if the reaction is an orphan, append its ID and index to the list
            orphan_reactions.append((reaction.id, index))
        index+=1

    # return the list of orphan reaction IDs
    return orphan_reactions


def find_exchange_rxns(model, inclObjFlag=False, irrevFlag=False):
    """
    find_exchange_rxns() find exchange and uptake reactions in a GEM model

    Args:
        - model (cobra.Model): The COBRA model.
        - inclObjFlag (bool, optional): Whether to include objective reactions in the exchange reaction set. Default to False.
        - irrevFlag (bool, optional): Whether the model is in irreversible format. Default to False.

    Returns:
        - tuple: boolean list indicating whether each reaction in the model is an exchange reaction (selExc),
                 and a boolean list indicating whether each reaction in the model is a nutrient uptake reaction (selUpt).
    """

    # create the stoichiometric matrix S from the model
    S = cobra.util.create_stoichiometric_matrix(model)

    if not irrevFlag:
        # find exchange reactions
        selExc = (S.sum(axis=0) == -1) | (S.sum(axis=0) == 1)
        selExc &= (S != 0).sum(axis=0) == 1

        if hasattr(model, 'objective_coefficients'):
            # remove objective reactions
            if not inclObjFlag:
                selExc[model.objective_coefficients != 0] = False
            else:
                selExc[model.objective_coefficients != 0] = True

        if hasattr(model, 'lower_bounds'):
            # find uptake reactions
            selUpt = model.lower_bounds < 0
            selUpt &= selExc
        else:
            selUpt = []

    else:
        # find exchange reactions
        selExc = (abs(S).sum(axis=0) == 1) & ((S != 0).sum(axis=0) == 1)

        if hasattr(model, 'objective_coefficients'):
            # remove objective reactions
            if not inclObjFlag:
                selExc[model.objective_coefficients != 0] = False
            else:
                selExc[model.objective_coefficients != 0] = True

        # find uptake reactions
        selUpt = (S.sum(axis=0) == 1) & ((S != 0).sum(axis=0) == 1)

    return selExc, selUpt


def add_reaction_annotations(input_file):
    data = []
    with open(input_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split('\t')
            rxn_id = line[0]
            rxn_name = line[1]
            annotation_dict = ast.literal_eval(line[2])
            data.append([rxn_id, rxn_name, annotation_dict])

    column_names = ['rxn_id', 'rxn_name']
    for d in data:
        for key in d[2].keys():
            if key not in column_names:
                column_names.append(key)

    df = pd.DataFrame(columns=column_names)
    for d in data:
        row_data = [d[0], d[1]]
        for col_name in column_names[2:]:
            if col_name in d[2]:
                row_data.append(d[2][col_name])
            else:
                row_data.append('')

        new_line = pd.DataFrame([row_data], columns=column_names)
        df = pd.concat([df, new_line])

    return df


def get_rxns_annotation(orphans, output):
    '''
    get_rxns_annotation() creates an output file containing the orphan reactions, their name and their annotation.
    '''

    results = []
    data = []

    # loop through the orphan reactions and write their information to the output file
    for rid in orphans:
        # get the reaction object from the model
        rxn = model.reactions.get_by_id(rid)
        
        # get the reaction name and annotation
        rxn_name = rxn.name
        rxn_annotation = rxn.annotation
        
        results.append((rid, rxn_name, str(rxn_annotation)))

    for r in results:
        rxn_id = r[0]
        rxn_name = r[1]
        annotation_dict = ast.literal_eval(r[2])
        data.append([rxn_id, rxn_name, annotation_dict])

    column_names = ['rxn_id', 'rxn_name']
    for d in data:
        for key in d[2].keys():
            if key not in column_names:
                column_names.append(key)

    df = pd.DataFrame(columns=column_names)
    for d in data:
        row_data = [d[0], d[1]]
        for col_name in column_names[2:]:
            if col_name in d[2]:
                row_data.append(d[2][col_name])
            else:
                row_data.append('')

        new_line = pd.DataFrame([row_data], columns=column_names)
        df = pd.concat([df, new_line])

    return df


# read in the SBML file as a COBRApy model
if args.extension == "xml" or args.extension == "sbml":
    model = cobra.io.read_sbml_model(args.input)
elif args.extension == "json":
    model = cobra.io.load_json_model(args.input)
elif args.extension == "mat":
    model = cobra.io.load_matlab_model(args.input)
else:
    raise ValueError("file extension not accepted")

# find all orphan reactions 
rxns = find_orphan_rxns(model)
indexes =  [i[1] for i in rxns]

# find the all the exchange in the model
selExc, selUpt = find_exchange_rxns(model, True, False)
non_exc_rxns = [idx for idx in indexes if selExc[idx-1]]

# find the orphan reactions (excluding exchange reactions)
orphans = set(indexes).difference(non_exc_rxns)
orphans_ids = [r[0] for r in rxns if r[1] in orphans]
print("Number of orphan reactions in the model: " + str(len(set(orphans))))

# get the annotations for the orphan reactions and write them to an output file
df = get_rxns_annotation(orphans_ids, args.output)
df.to_csv(args.output, sep='\t', index=False)


