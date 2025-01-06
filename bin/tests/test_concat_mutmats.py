import pytest

import pandas as pd 
from tempfile import TemporaryDirectory

import os
import sys
current_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(parent_dir)

from concat_mutmats import get_combined_mutmat

def test_get_combined_mutmat():
    indir = 'tests/test_data/mutation_matrices'
    # result = load_mutmat(indir=indir)
    # result_to_set = {k:set(result[k]) for k in result}
    # print(result)
    expected = pd.DataFrame({
        'MutationType': ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'],
        'PD47151_1': [261, 85, 888, 136, 274, 97],
        'PD47151_2': [242, 89, 993, 136, 259, 112],
        'PD47151_3': [209, 89, 745, 115, 208, 240],
        'PD52103_1': [263, 103, 739, 160, 305, 196],
        'PD52103_10': [633, 171, 1388, 226, 440, 317]
    })
    with TemporaryDirectory() as outdir:
        get_combined_mutmat(indir, outdir)
        result = pd.read_csv(f'{outdir}/combined_mutmat.SBS6.all', sep = '\t')
        print(result)
        result = result.reindex(columns = expected.columns)
        assert result.equals(expected)
    