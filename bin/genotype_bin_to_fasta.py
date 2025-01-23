#!/usr/bin/env python3

import pandas as pd

import argparse, sys

def genotype_bin_to_seq(binary, ref, alt):
    switcher = {
        0: ref,
        1: alt,
        0.5: '?'
    }
    return switcher[binary]

def genotype_bin_to_fasta(genotype_bin_path, outfile):
    genotype_bin = pd.read_csv(genotype_bin_path, sep = ' ')
    (_, _, ref, alt) = list(zip(*(m.split('_') for m in genotype_bin.index)))
    with open(outfile, mode='w') as fasta:
        refseq = ''.join(ref)
        fasta.write(f'>Ancestral\n{refseq}\n')
        for sample in genotype_bin.columns:
            seq = ''.join(list(map(genotype_bin_to_seq, genotype_bin[sample], ref, alt)))
            fasta.write(f'>{sample}\n{seq}\n')
            

def get_arguments():
    parser = argparse.ArgumentParser(description='Binary genotype to fasta')
    parser.add_argument('--genotype_bin_path', required=True,
                        help='path to genotype bin file', type = str)
    parser.add_argument('--outfile', '-o', required=True,
                        help='output file', type = str)
    return parser

def main(args):
    genotype_bin_to_fasta(args.genotype_bin_path, args.outfile)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))