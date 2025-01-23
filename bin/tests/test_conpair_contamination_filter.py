import pytest

import pandas as pd 

import os
import sys
current_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(parent_dir)

from conpair_contamination_filter import filter_by_concordance, get_contaminated_samples



def test_filter_by_concordance():
    samples_dict = {'PD47151n_lo0002': ['PD47151n_lo0002',
        'PD47151b',
        'PD47151',
        'some/path/47151'],
        'PD47151n_lo0004': ['PD47151n_lo0004',
        'PD47151b',
        'PD47151',
        'some/path/47151_4'],
        'PD52103n_lo0002': ['PD52103n_lo0002',
        'PD52103b',
        'PD52103',
        'some/path/52103']}
    concordance = pd.DataFrame([['PD52103n_lo0002', 'PD47151b', 24.65, 0.5566435468516252],
       ['PD47151n_lo0004', 'PD52103b', 25.72, 0.4736842105263157],
       ['PD52103n_lo0002', 'PD52103b', 99.8, 0.5419556643546851],
       ['PD47151n_lo0002', 'PD47151b', 99.61, 0.4547803617571059],
       ['PD47151n_lo0002', 'PD52103b', 24.92, 0.4432204542363661],
       ['PD47151n_lo0004', 'PD47151b', 99.78, 0.489188086495308]],
        columns=['sample_id', 'match_id', 'concordance', 'fraction_of_markers'])
    
    expected = {'PD52103n_lo0002': ['PD52103b'],
        'PD47151n_lo0002': ['PD47151b'],
        'PD47151n_lo0004': ['PD47151b']}
    assert filter_by_concordance(samples_dict, concordance, concordance_threshold = 90) == expected

def test_filter_by_concordance_one_wrong_match(caplog):
    samples_dict = {'PD47151n_lo0002': ['PD47151n_lo0002',
        'PD47151b',
        'PD47151',
        'some/path/47151'],
        'PD47151n_lo0004': ['PD47151n_lo0004',
        'PD47151b',
        'PD47151',
        'some/path/47151_4'],
        'PD52103n_lo0002': ['PD52103n_lo0002',
        'PD52103b',
        'PD52103',
        'some/path/52103']}
    concordance = pd.DataFrame([['PD52103n_lo0002', 'PD47151b', 24.65, 0.5566435468516252],
       ['PD47151n_lo0004', 'PD52103b', 95, 0.4736842105263157],
       ['PD52103n_lo0002', 'PD52103b', 99.8, 0.5419556643546851],
       ['PD47151n_lo0002', 'PD47151b', 99.61, 0.4547803617571059],
       ['PD47151n_lo0002', 'PD52103b', 24.92, 0.4432204542363661],
       ['PD47151n_lo0004', 'PD47151b', 20, 0.489188086495308]],
        columns=['sample_id', 'match_id', 'concordance', 'fraction_of_markers'])
    
    expected = {'PD52103n_lo0002': ['PD52103b'],
        'PD47151n_lo0002': ['PD47151b']}
    assert filter_by_concordance(samples_dict, concordance, concordance_threshold = 90) == expected
    assert "sample PD47151n_lo0004 matches the wrong match normal ['PD52103b'], potentially a labelling error,\nplease check whether this is the case,\nremoving sample PD47151n_lo0004" in caplog.text

def test_filter_by_concordance_one_multiple_match(caplog):
    samples_dict = {'PD47151n_lo0002': ['PD47151n_lo0002',
        'PD47151b',
        'PD47151',
        'some/path/47151'],
        'PD47151n_lo0004': ['PD47151n_lo0004',
        'PD47151b',
        'PD47151',
        'some/path/47151_4'],
        'PD52103n_lo0002': ['PD52103n_lo0002',
        'PD52103b',
        'PD52103',
        'some/path/52103']}
    concordance = pd.DataFrame([['PD52103n_lo0002', 'PD47151b', 24.65, 0.5566435468516252],
       ['PD47151n_lo0004', 'PD52103b', 95, 0.4736842105263157],
       ['PD52103n_lo0002', 'PD52103b', 99.8, 0.5419556643546851],
       ['PD47151n_lo0002', 'PD47151b', 99.61, 0.4547803617571059],
       ['PD47151n_lo0002', 'PD52103b', 24.92, 0.4432204542363661],
       ['PD47151n_lo0004', 'PD47151b', 99.78, 0.489188086495308]],
        columns=['sample_id', 'match_id', 'concordance', 'fraction_of_markers'])
    
    expected = {'PD52103n_lo0002': ['PD52103b'],
        'PD47151n_lo0002': ['PD47151b']}
    assert filter_by_concordance(samples_dict, concordance, concordance_threshold = 90) == expected
    assert "sample PD47151n_lo0004 matches more than one match normal ['PD52103b', 'PD47151b'], samples were potentially mixed up, \nplease check whether this is the case, \nremoving sample PD47151n_lo0004" in caplog.text
    
def test_filter_by_concordance_one_no_match(caplog):
    samples_dict = {'PD47151n_lo0002': ['PD47151n_lo0002',
        'PD47151b',
        'PD47151',
        'some/path/47151'],
        'PD47151n_lo0004': ['PD47151n_lo0004',
        'PD47151b',
        'PD47151',
        'some/path/47151_4'],
        'PD52103n_lo0002': ['PD52103n_lo0002',
        'PD52103b',
        'PD52103',
        'some/path/52103']}
    concordance = pd.DataFrame([['PD52103n_lo0002', 'PD47151b', 24.65, 0.5566435468516252],
       ['PD47151n_lo0004', 'PD52103b', 25.72, 0.4736842105263157],
       ['PD52103n_lo0002', 'PD52103b', 99.8, 0.5419556643546851],
       ['PD47151n_lo0002', 'PD47151b', 99.61, 0.4547803617571059],
       ['PD47151n_lo0002', 'PD52103b', 24.92, 0.4432204542363661],
       ['PD47151n_lo0004', 'PD47151b', 30, 0.489188086495308]],
        columns=['sample_id', 'match_id', 'concordance', 'fraction_of_markers'])
    
    expected = {'PD52103n_lo0002': ['PD52103b'],
        'PD47151n_lo0002': ['PD47151b']}
    assert filter_by_concordance(samples_dict, concordance, concordance_threshold = 90) == expected
    assert "samples ['PD47151n_lo0004'] did not match any match normal" in caplog.text
    
def test_get_contaminated_samples():
    contamination = pd.DataFrame(
        [['PD52103n_lo0002', 'PD52103b', 0.094, 0.015],
         ['PD47151n_lo0002', 'PD47151b', 0.5, 0.025]],
        columns = ['sample_id', 'match_id', 'contamination_sample', 'contamination_match']
    )
    expected = {"contaminated_samples": ['PD47151n_lo0002'], "contaminated_matches": []}
    assert get_contaminated_samples(contamination) == expected
    
def test_get_contaminated_match():
    contamination = pd.DataFrame(
        [['PD52103n_lo0002', 'PD52103b', 0.094, 10],
         ['PD47151n_lo0002', 'PD47151b', 0.5, 0.025]],
        columns = ['sample_id', 'match_id', 'contamination_sample', 'contamination_match']
    )
    expected = {"contaminated_samples": ['PD47151n_lo0002'], "contaminated_matches": ['PD52103b']}
    assert get_contaminated_samples(contamination) == expected