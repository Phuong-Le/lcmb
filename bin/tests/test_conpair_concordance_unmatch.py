import pytest

import pandas as pd

import os
import sys
current_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(parent_dir)

from conpair_concordance_filter_unmatch import get_problematic_concordance



def test_get_problematic_concordance_pass(caplog):
    samples_pdid_dict = {'PD47151n_lo0002': 'PD47151',
        'PD47151n_lo0004': 'PD47151',
        'PD52103n_lo0002': 'PD52103' }
    concordance = pd.DataFrame([['PD52103n_lo0002', 'PD47151n_lo0002', 24.65, 0.5566435468516252],
    ['PD52103n_lo0002', 'PD47151n_lo0004', 25, 0.4736842105263157],
    ['PD47151n_lo0004', 'PD47151n_lo0002', 91, 0.4736842105263157],
    ['PD47151n_lo0004', 'PD52103n_lo0002', 25, 0.4736842105263157],
    ['PD47151n_lo0002', 'PD47151n_lo0004', 99.61, 0.4547803617571059],
    ['PD47151n_lo0002', 'PD52103n_lo0002', 24.92, 0.4432204542363661]],
    columns=['sample_id', 'match_id', 'concordance', 'fraction_of_markers'])

    expected = {}
    assert get_problematic_concordance(concordance, samples_pdid_dict, 90) == expected
    assert "cannot check {'PD52103n_lo0002'} because they do not match any other samples" in caplog.text


def test_get_problematic_concordance_wrong_match(caplog):
    samples_pdid_dict = {'PD47151n_lo0002': 'PD47151',
        'PD47151n_lo0004': 'PD47151',
        'PD52103n_lo0002': 'PD52103' }
    concordance = pd.DataFrame([['PD52103n_lo0002', 'PD47151n_lo0002', 92, 0.5566435468516252],
    ['PD52103n_lo0002', 'PD47151n_lo0004', 25, 0.4736842105263157],
    ['PD47151n_lo0004', 'PD47151n_lo0002', 22, 0.4736842105263157],
    ['PD47151n_lo0004', 'PD52103n_lo0002', 25, 0.4736842105263157],
    ['PD47151n_lo0002', 'PD47151n_lo0004', 25, 0.4547803617571059],
    ['PD47151n_lo0002', 'PD52103n_lo0002', 92, 0.4432204542363661]],
    columns=['sample_id', 'match_id', 'concordance', 'fraction_of_markers'])

    expected = {'PD52103n_lo0002': ['PD47151n_lo0002'], 'PD47151n_lo0002': ['PD52103n_lo0002']}
    assert get_problematic_concordance(concordance, samples_pdid_dict, 90) == expected
    assert "sample PD52103n_lo0002 matches the wrong match normal ['PD47151n_lo0002'], \nremoving sample PD52103n_lo0002" in caplog.text
    assert "sample PD47151n_lo0002 matches the wrong match normal ['PD52103n_lo0002'], \nremoving sample PD47151n_lo0002" in caplog.text


def test_get_problematic_concordance_multiple_match(caplog):
    samples_pdid_dict = {'PD47151n_lo0002': 'PD47151',
        'PD47151n_lo0004': 'PD47151',
        'PD52103n_lo0002': 'PD52103' }
    concordance = pd.DataFrame([['PD52103n_lo0002', 'PD47151n_lo0002', 92, 0.5566435468516252],
    ['PD52103n_lo0002', 'PD47151n_lo0004', 25, 0.4736842105263157],
    ['PD47151n_lo0004', 'PD47151n_lo0002', 92, 0.4736842105263157],
    ['PD47151n_lo0004', 'PD52103n_lo0002', 25, 0.4736842105263157],
    ['PD47151n_lo0002', 'PD47151n_lo0004', 91, 0.4547803617571059],
    ['PD47151n_lo0002', 'PD52103n_lo0002', 92, 0.4432204542363661]],
    columns=['sample_id', 'match_id', 'concordance', 'fraction_of_markers'])

    expected = {'PD47151n_lo0002': ['PD47151n_lo0004', 'PD52103n_lo0002'],
                'PD52103n_lo0002': ['PD47151n_lo0002']}
    assert get_problematic_concordance(concordance, samples_pdid_dict, 90) == expected
    assert "sample PD47151n_lo0002 matches more than one match normal ['PD47151n_lo0004', 'PD52103n_lo0002'], \nremoving sample PD47151n_lo0002" in caplog.text
