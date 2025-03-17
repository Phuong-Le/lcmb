#!/usr/bin/env python3

from pysam import tabix_compress, tabix_index

import argparse, sys

def get_arguments():
    parser = argparse.ArgumentParser(description='tabix compress and index a file')
    parser.add_argument('--infile', required=True,
                        help='input file', type = str)
    parser.add_argument('--output_gz_path', '-o', required=True,
                        help='output bgzipped file', type = str)
    parser.add_argument('--preset', required=True,
                        help='preset, Valid values for preset are gff, bed, sam, vcf, psltbl, pileup', type = str)
    return parser

def main(args):
    tabix_compress(args.infile, args.output_gz_path)
    tabix_index(args.output_gz_path, preset=args.preset)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
