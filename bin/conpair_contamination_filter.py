#!/usr/bin/env python3

import pandas as pd
from collections import defaultdict
import logging


import argparse, sys

logging.basicConfig(filename='conpair_filter.log', filemode='w', level = logging.INFO)

def filter_by_concordance(samples_dict, concordance, concordance_threshold = 90):
    """take a concordance dictionary and filter out samples that match to exactly one match normal

    Args:
        samples_dict (dict): sample dictionary with key as sample ID and values being all elements in the sample paths file
        concordance (dataframe): dataframe that contains concordance information between samples and all match normal
        concordance_threshold (float, optional): concordance threshold taking values from 0-100. Defaults to 90.

    Returns:
        dict: a filtered dictionary for concordance
    """
    concordance_filtered = concordance[concordance['concordance'] > concordance_threshold]
    concordance_filtered_list = concordance_filtered[['sample_id', 'match_id']].to_numpy().tolist()
    concordance_dict = defaultdict(list)
    for sample, match in concordance_filtered_list:
        concordance_dict[sample].append(match)
    concordance_dict_unique = {}
    for k,v in concordance_dict.items():
        if len(v) == 1:
            if concordance_dict[k][0] == samples_dict[k][1]:
                concordance_dict_unique[k] = v
            else:
                logging.warning(f'sample {k} matches the wrong match normal {v}, potentially a labelling error,\nplease check whether this is the case,\nremoving sample {k}')
        else:
            logging.warning(f'sample {k} matches more than one match normal {v}, samples were potentially mixed up, \nplease check whether this is the case, \nremoving sample {k}')
    no_match = [sample for sample in samples_dict if sample not in concordance_dict]
    if no_match:
        logging.warning(f'samples {no_match} did not match any match normal')
    return concordance_dict_unique


def get_contaminated_samples(contamination, contamination_threshold_samples = 0.3, contamination_threshold_match = 5):
    """identify contaminated samples based on the threshold

    Args:
        contamination (dataframe): contamination dataframe with columns being ['sample_id', 'match_id', 'contamination_sample', 'contamination_match'], the last column indicates contamination score
        contamination_threshold (float, optional): contamination threhold. Defaults to 0.1.

    Returns:
        dict: contaminated samples for the samples and the match normal
    """
    logging.info(f'contamination threhold for the samples is {contamination_threshold_samples}')
    logging.info(f'contamination threhold for the match normal is {contamination_threshold_match}')
    # contamination for the samples
    sample_contamination = dict(contamination[['sample_id', 'contamination_sample']].to_numpy().tolist())
    contaminated_samples = [sample for sample in sample_contamination if sample_contamination[sample] > contamination_threshold_samples]
    # contamination for the match normal
    match_contamination = dict(contamination[['match_id', 'contamination_match']].to_numpy().tolist())
    contaminated_matches = [sample for sample in match_contamination if match_contamination[sample] > contamination_threshold_match]
    return {"contaminated_samples": contaminated_samples, "contaminated_matches": contaminated_matches}

def filter_contaminations(samples_path, concordance_path, contamination_path, outfile, concordance_threshold = 90, contamination_threshold_samples = 0.3, contamination_threshold_match = 5):
    samples = pd.read_csv(samples_path, sep = '\t')
    samples_dict = {row[0]: row for row in samples.to_numpy().tolist()}

    concordance = pd.read_csv(concordance_path, sep = '\t', names = ['sample_id', 'match_id', 'concordance', 'fraction_of_markers'])
    contamination = pd.read_csv(contamination_path, sep = '\t', names = ['sample_id', 'match_id', 'contamination_sample', 'contamination_match'])

    # concordance
    concordance_dict_unique = filter_by_concordance(samples_dict, concordance, concordance_threshold)
    samples_concordance_filtered = samples[samples['sample_id'].isin(concordance_dict_unique.keys())]
    # contamination for the samples
    contaminated_samples = get_contaminated_samples(contamination, contamination_threshold_samples = contamination_threshold_samples, contamination_threshold_match = contamination_threshold_match)
    logging.warning(f'removing {contaminated_samples} as they are contaminated')
    samples_concordance_contamination_filtered = samples_concordance_filtered[(~samples_concordance_filtered['sample_id'].isin(contaminated_samples["contaminated_samples"])) & (~samples_concordance_filtered['match_normal_id'].isin(contaminated_samples["contaminated_matches"]))]

    samples_concordance_contamination_filtered.to_csv(outfile, index = False, sep = '\t')
    logging.info(f'the filtered sample paths are now stored in {outfile}')

def get_arguments():
    parser = argparse.ArgumentParser(description='Check concordance and contamination')
    parser.add_argument('--samples_path', required=True,
                        help='path to the file that contains paths to the samples', type = str)
    parser.add_argument('--concordance_path', required=True,
                        help='path to the concordance file', type = str)
    parser.add_argument('--contamination_path', required=True,
                        help='path to the contamination file', type = str)
    parser.add_argument('--outfile', '-o', required=True,
                        help='output file', type = str)
    parser.add_argument('--concordance_threshold', required=False,
                        help='concordance threshold, default to 90', type = float, default = 90)
    parser.add_argument('--contamination_threshold_samples', required=False,
                        help='contamination threshold for the samples', type = float, default = 0.3)
    parser.add_argument('--contamination_threshold_match', required=False,
                        help='contamination threshold for the match normal', type = float, default = 5)
    return parser

def main(args):
    filter_contaminations(args.samples_path, args.concordance_path, args.contamination_path, args.outfile, concordance_threshold = args.concordance_threshold, contamination_threshold_samples = args.contamination_threshold_samples, contamination_threshold_match = args.contamination_threshold_match)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))


