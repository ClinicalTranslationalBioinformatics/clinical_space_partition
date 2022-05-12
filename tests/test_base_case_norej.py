import pytest

from csp import csp_norej, obtain_predictor_intervals


@pytest.fixture
def base_case():
    return {
        'filename': '../demo/csp-norej.config',
        'rho': 0.5,
        'predictors': {
            'CADD': [0.995, 0.254],
            'MutPred': [0.95, 0.706],
            'PolyPhen-2': [0.926, 0.638],
            'SIFT': [0.924, 0.682],
            'VEST': [0.971, 0.824],
            'fathmm': [0.835, 0.658]
        },
        'x_points': [
            0.040404040404040435,
            0.043478260869565244,
            0.09054325955734414,
            0.14228456913827642,
            0.15231788079470185,
            0.28368794326241126,
            0.819819819819819
        ],
        'interval_best_predictors': {
            (0, 0.040404040404040435): 'CADD',
            (0.040404040404040435, 0.043478260869565244): 'VEST',
            (0.043478260869565244, 0.09054325955734414): 'VEST',
            (0.09054325955734414, 0.14228456913827642): 'VEST',
            (0.14228456913827642, 0.15231788079470185): 'VEST',
            (0.15231788079470185, 0.28368794326241126): 'VEST',
            (0.28368794326241126, 0.819819819819819): 'VEST',
            (0.819819819819819, 1): 'VEST'
        },
        'merged_intervals': {
            'CADD': 0.040404040404040435,
            'VEST': 0.9595959595959596
        },
        'output': '''
CLINICAL SPACE PARTITION
------------------------

Parameters considered: sensitivity and specificity

Methods compared: CADD, MutPred, PolyPhen-2, SIFT, VEST, fathmm

Best combination of methods (rho=0.5): VEST, CADD

List of clinical space fraction for each predictor:

Predictor 	Relative value
--------- 	--------------
VEST      	0.96
CADD      	0.04
MutPred   	0
PolyPhen-2	0
SIFT      	0
fathmm    	0
'''
    }


def test_parse_config(base_case):
    rho, predictors = csp_norej.parse_config(base_case['filename'], mode='norej')
    assert rho == base_case['rho']
    assert predictors == base_case['predictors']


def test_main(base_case, capsys):
    csp_norej.main(base_case['rho'], base_case['predictors'])
    captured = capsys.readouterr()
    assert captured.out == base_case['output']


def test_get_predictors_intersection(base_case):
    x_points = obtain_predictor_intervals.get_predictors_intersections(base_case['rho'], base_case['predictors'])
    assert x_points == base_case['x_points']


def test_get_interval_best_predictor(base_case):
    interval_best_predictors = obtain_predictor_intervals.get_interval_best_predictor(base_case['rho'],
                                                                                      base_case['predictors'],
                                                                                      base_case['x_points'])
    assert interval_best_predictors == base_case['interval_best_predictors']


def test_merge_intervals(base_case):
    merged_intervals = obtain_predictor_intervals.merge_intervals(base_case['interval_best_predictors'])
    assert merged_intervals == base_case['merged_intervals']
