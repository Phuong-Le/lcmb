#!/usr/bin/env python3

import pandas as pd
from collections import defaultdict
from pathlib import Path
from glob import glob

import argparse, sys


def load_mutmat(indir):
    """recursively search a directory for files with extension ".all", 
    these are assumed to be mutation matrices

    Args:
        indir (path): directory to be searched

    Returns:
        dict: a dictionary with key as file extensions, and values as a list of dataframes from files with these extensions
    """
    mutmat_paths = glob(f'{indir}/**/*.all', recursive=True)
    mutmat_dict = defaultdict(list)
    for mutmat_path in mutmat_paths:
        mutmat = pd.read_csv(mutmat_path, sep = '\t')
        context_type = ''.join(Path(mutmat_path).suffixes) # context type based on extensions 
        mutmat_dict[context_type].append(mutmat)
    return mutmat_dict


def concat_mutmat(mutmats):
    """concatenating a list of mutation matrix into one mutation matrix

    Args:
        mutmats (list): list of mutation matrices

    Returns:
        DataFrame: a dataframe of the combined mutation matrix
    """
    mut_vector = mutmats[0]['MutationType']
    corrected_mutmats = [mut_vector]
    for mutmat in mutmats:
        mutmat = mutmat.set_index('MutationType', drop=True) 
        mutmat = mutmat.reindex(mut_vector)
        mutmat.reset_index(drop=True, inplace=True)
        corrected_mutmats.append(mutmat)
    return pd.concat(corrected_mutmats, axis=1)


def get_combined_mutmat(indir, outdir):
    """output combined mutation matrices given a directory that contains mutation matrices, 
    assuming the matrices to be combined have the same extensions

    Args:
        indir (path): path to the directory where the mutation matrix with extensions '.all' will be recursively searched for
        outdir (path): path to the directory where the combined mutation matrix of all files with the same extensions are written to
    """
    mutmat_dict = load_mutmat(indir)
    # print(mutmat_dict)
    for context_type in mutmat_dict:
        combined_mutmat = concat_mutmat(mutmat_dict[context_type])
        # print(combined_mutmat)
        outpath = f'{outdir}/combined_mutmat{context_type}'
        combined_mutmat.to_csv(outpath, sep = '\t', index = False)


def get_arguments():
    parser = argparse.ArgumentParser(description='Combining mutation matrices for the branches across donors')
    parser.add_argument('--indir', required=True,
                        help="path to the directory where the mutation matrix with extensions '.all' will be recursively searched for", type = str)
    parser.add_argument('--outdir', required=True,
                        help="path to the directory where the combined mutation matrix of all files with the same extensions are written to", type = str)
    return parser


def main(args):
    get_combined_mutmat(args.indir, args.outdir)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))

            
