import pytest

import pandas as pd 
from tempfile import NamedTemporaryFile

import os
import sys
current_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(parent_dir)

from genotype_bin_to_fasta import genotype_bin_to_fasta

def test_genotype_bin_to_fasta():
    genotype_bin = pd.DataFrame(
         [[0, 1, 1],
        [0, 0, 1],
        [1, 0, 1]],
         columns= ['sample1', 'sample2', 'sample3'],
         index = ['chr1_102_G_A', 'chr4_105_T_C', 'chr10_102_G_C']
    )
    with NamedTemporaryFile(delete=False, mode = 'w+t') as genotype_bin_path:    
          genotype_bin.to_csv(genotype_bin_path.name, sep = ' ', index_label = False)
          with NamedTemporaryFile(delete=False, mode = 'w+t') as fasta:
               genotype_bin_to_fasta(genotype_bin_path, fasta.name)
               expected = ">Ancestral\nGTG\n>sample1\nGTC\n>sample2\nATG\n>sample3\nACC\n"
               with open(fasta.name) as outfasta:
                    content = outfasta.read()
                    print(content == expected)