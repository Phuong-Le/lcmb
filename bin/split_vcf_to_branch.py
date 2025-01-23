#!/usr/bin/env python3

import pandas as pd
import subprocess
import os

import argparse, sys


def split_vcf(vcf_path, outdir, prefix):
    vcf = pd.read_csv(vcf_path, sep = '\t')
    vcf = vcf[['CHROM', 'POS', 'Branch', 'REF', 'ALT']] # rearrange because column index is what MatrixGenerator uses to generate mutation matrix
    branches = vcf['Branch'].unique()
    os.mkdir(f'{outdir}/vcf_with_header')
    os.mkdir(f'{outdir}/matrix_generator')
    for branch in branches:
        vcf_branch = vcf[vcf['Branch'] == branch]
        outfile_with_header = f'{outdir}/vcf_with_header/{prefix}_{branch}_with_header.txt'
        outfile = f'{outdir}/matrix_generator/{prefix}_{branch}.vcf'
        vcf_branch.to_csv(outfile_with_header, index = False, sep = '\t')
        cmd = f"sed '1s/^/#/' {outfile_with_header} > {outfile}"
        subprocess.run(cmd, shell=True)
    
def get_arguments():
    parser = argparse.ArgumentParser(description='Binary genotype to fasta')
    parser.add_argument('--vcf_path', required=True,
                        help='path to big vcf', type = str)
    parser.add_argument('--outdir', '-o', required=True,
                        help='output directory', type = str)
    parser.add_argument('--prefix', required=True,
                        help='output prefix', type = str)
    return parser

def main(args):
    split_vcf(args.vcf_path, args.outdir, args.prefix)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
