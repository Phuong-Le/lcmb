#!/usr/bin/env python3

import pandas as pd
from collections import defaultdict
import logging


import argparse, sys

logging.basicConfig(filename='conpair_filter.log', filemode='w', level = logging.INFO)


def filter_contamination_unmatch(samples_path, contamination_path, outfile, contamination_threshold_samples):
    samples = pd.read_csv(samples_path, sep = '\t')
    contamination = pd.read_csv(contamination_path, sep = '\t', names = ['sample_id', 'match_id', 'contamination_sample', 'contamination_match'])
    available_samples = contamination['sample_id']
    contaminated_samples = contamination[contamination['contamination_sample'] > contamination_threshold_samples]

    unchecked_samples = [sample for sample in samples['sample_id'].tolist() if sample not in available_samples.tolist()]
    logging.warn(f"cannot check {unchecked_samples} because there are no other samples to compare to")

    samples_contamination_filtered = samples[(~samples['sample_id'].isin(contaminated_samples))]
    samples_contamination_filtered.to_csv(outfile, index = False, sep = '\t')
    logging.info(f'the filtered sample paths are now stored in {outfile}')


def get_arguments():
    parser = argparse.ArgumentParser(description='Check contamination')
    parser.add_argument('--samples_path', required=True,
                        help='path to the file that contains paths to the samples', type = str)
    parser.add_argument('--contamination_path', required=True,
                        help='path to the contamination file', type = str)
    parser.add_argument('--outfile', '-o', required=True,
                        help='output file', type = str)
    parser.add_argument('--contamination_threshold_samples', required=False,
                        help='contamination threshold for the samples', type = float, default = 0.3)

    return parser

def main(args):
    filter_contamination_unmatch(args.samples_path, args.contamination_path, args.outfile, args.contamination_threshold_samples)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))


