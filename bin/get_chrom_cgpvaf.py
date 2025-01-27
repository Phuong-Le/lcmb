#!/usr/bin/env python3

import pandas as pd

import argparse, sys


def get_chrom_list(fai, hdr, outfile):
    fai = pd.read_csv(fai, sep="\t", header=None, names = ["chrom", "length", "offset", "linebases", "linewidth"])
    hdr = pd.read_csv(hdr, sep="\t", header=None, usecols=[0,1,2], names=["chrom", "start", "end"])
    hdr_chrom_to_exclude = hdr[hdr[['chrom', 'end']].apply(lambda x: x.to_numpy().tolist() in fai[['chrom', 'length']].to_numpy().tolist(), axis=1)]
    chrom_list = fai[~fai['chrom'].isin(hdr_chrom_to_exclude['chrom'])]['chrom']
    chrom_list.to_csv(outfile, sep="\t", header=False, index=False)

def get_arguments():
    parser = argparse.ArgumentParser(description='Get relevant chromosome list for cgpVAF based on fasta index and excluding high depth regions')
    parser.add_argument('--fai', required=True,
                        help="path to fai", type = str)
    parser.add_argument('--hdr', required=True,
                        help="path to high depth region", type = str)
    parser.add_argument('--outfile', required=True,
                        help="path to output chrom list file", type = str)
    return parser


def main(args):
    get_chrom_list(args.fai, args.hdr, args.outfile)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))


