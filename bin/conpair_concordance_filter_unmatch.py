#!/usr/bin/env python3

import pandas as pd
from collections import defaultdict
import logging


import argparse, sys

logging.basicConfig(filename='conpair_filter.log', filemode='w', level = logging.INFO)

def get_problematic_concordance(concordance, samples_pdid_dict, concordance_threshold):
    concordance_filtered = concordance[concordance['concordance'] > concordance_threshold]
    concordance_filtered_list = concordance_filtered[['sample_id', 'match_id']].to_numpy().tolist()

    concordance_dict = defaultdict(list)
    concordance_pdid_dict = defaultdict(list)
    for sample, match in concordance_filtered_list:
        concordance_dict[sample].append(match)
        concordance_pdid_dict[sample].append(samples_pdid_dict[match])

    no_match = {sample_id for sample_id in concordance['sample_id'] if sample_id not in concordance_dict}
    logging.warning(f'cannot check {no_match} because they do not match any other samples')

    concordance_problematic = {}
    for k, v in concordance_dict.items():
        match_pdids = list(set(concordance_pdid_dict[k]))
        if len(match_pdids) == 1:
            if samples_pdid_dict[k] != match_pdids[0]:
                logging.warning(f'sample {k} matches the wrong match normal {v}, \nremoving sample {k}')
                concordance_problematic[k] = v
        elif len(match_pdids) > 1:
            logging.warning(f'sample {k} matches more than one match normal {v}, \nremoving sample {k}')
            concordance_problematic[k] = v

    return concordance_problematic

def filter_concordance_unmatch(samples_path, concordance_path, outfile, concordance_threshold = 90):
    samples = pd.read_csv(samples_path, sep = '\t')
    sample_ids = samples['sample_id']
    pdids = samples['pdid']
    samples_pdid_dict = {sample_ids[i]: pdids[i] for i in range(len(sample_ids))}

    concordance = pd.read_csv(concordance_path, sep = '\t', names = ['sample_id', 'match_id', 'concordance', 'fraction_of_markers'])

    concordance_problematic = get_problematic_concordance(concordance, samples_pdid_dict, concordance_threshold)
    samples_concordance_filtered = samples[(~samples['sample_id'].isin(concordance_problematic.keys()))]
    samples_concordance_filtered.to_csv(outfile, index = False, sep = '\t')
    logging.info(f'the filtered sample paths are now stored in {outfile}')
    return samples_concordance_filtered

def get_arguments():
    parser = argparse.ArgumentParser(description='Check concordance and contamination')
    parser.add_argument('--samples_path', required=True,
                        help='path to the file that contains paths to the samples', type = str)
    parser.add_argument('--concordance_path', required=True,
                        help='path to the concordance file', type = str)
    parser.add_argument('--outfile', '-o', required=True,
                        help='output file', type = str)
    parser.add_argument('--concordance_threshold', required=False,
                        help='concordance threshold, default to 90', type = float, default = 90)
    return parser

def main(args):
    filter_concordance_unmatch(args.samples_path, args.concordance_path, args.outfile, concordance_threshold = args.concordance_threshold)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))


