import pytest

from csp import csp_co, find_predictor_intersections, build_intersection_graph, search_graph_polygons, \
    obtain_polygon_data


@pytest.fixture
def base_case():
    return {
        'filename': '../demo/csp-co.config',
        'rho': 0.5,
        'predictors_demo': {
            'CADD': [0.995, 0.254, 1.0],
            'MutPred': [0.95, 0.706, 0.281],
            'PolyPhen-2': [0.926, 0.638, 0.909],
            'SIFT': [0.924, 0.682, 0.866],
            'VEST': [0.971, 0.824, 0.937],
            'fathmm': [0.835, 0.658, 0.902]
        },
        'predictors': {
            'PolyPhen-2': [
                0.926,
                0.638,
                0.909
            ],
            'SIFT': [
                0.924,
                0.682,
                0.866
            ],
            'CADD': [
                0.995,
                0.254,
                1
            ]
        },
        'nodes': {
            (0.0, 0.0): {
                'x_axis',
                'y_axis'
            },
            (0.0, 1.0): {
                'y_axis',
                'hypotenuse'
            },
            (1.0, 0.0): {
                'hypotenuse',
                'x_axis'
            },
            (0.98341909, 0.0): {
                'x_axis',
                (-0.62611871, 0.61573709)
            },
            (0.0, 0.61573709): {
                'y_axis',
                (-0.62611871, 0.61573709)
            },
            (0.0, 0.30386916): {
                'y_axis',
                (-0.19990917, 0.30386916)
            },
            (0.0, 0.36284274): {
                'y_axis',
                (-0.28050451, 0.36284274)
            },
            (0.87006477, 0.12993523): {
                'hypotenuse',
                (-0.19990917, 0.30386916)
            },
            (0.88556116, 0.11443884): {
                'hypotenuse',
                (-0.28050451, 0.36284274)
            },
            (0.731724419887841, 0.157590740144327): {
                (-0.62611871, 0.61573709),
                (-0.28050451, 0.36284274)
            },
            (0.731724423625055, 0.157590737804387): {
                (-0.62611871, 0.61573709),
                (-0.19990917, 0.30386916)
            },
            (0.731724439651225, 0.157590734600608): {
                (-0.28050451, 0.36284274),
                (-0.19990917, 0.30386916)
            }
        },
        'lines': {
            'x_axis': [
                (0.0, 0.0),
                (0.98341909, 0.0),
                (1.0, 0.0)
            ],
            'y_axis': [
                (0.0, 0.0),
                (0.0, 0.30386916),
                (0.0, 0.36284274),
                (0.0, 0.61573709),
                (0.0, 1.0)
            ],
            'hypotenuse': [
                (0.0, 1.0),
                (0.87006477, 0.12993523),
                (0.88556116, 0.11443884),
                (1.0, 0.0)],
            (-0.62611871, 0.61573709): [
                (0.0, 0.61573709),
                (0.731724419887841, 0.157590740144327),
                (0.731724423625055, 0.157590737804387),
                (0.98341909, 0.0)],
            (-0.19990917, 0.30386916): [
                (0.0, 0.30386916),
                (0.731724423625055, 0.157590737804387),
                (0.731724439651225, 0.157590734600608),
                (0.87006477, 0.12993523)
            ],
            (-0.28050451, 0.36284274): [
                (0.0, 0.36284274),
                (0.731724419887841, 0.157590740144327),
                (0.731724439651225, 0.157590734600608),
                (0.88556116, 0.11443884)
            ]
        },
        'interactions': {
            (0.0, 0.0): [
                (0.0, 0.30386916),
                (0.98341909, 0.0)
            ],
            (0.0, 0.30386916): [
                (0.0, 0.0),
                (0.0, 0.36284274),
                (0.731724423625055, 0.157590737804387)
            ],
            (0.0, 0.36284274): [
                (0.0, 0.30386916),
                (0.0, 0.61573709),
                (0.731724419887841, 0.157590740144327)
            ],
            (0.0, 0.61573709): [
                (0.0, 0.36284274),
                (0.0, 1.0),
                (0.731724419887841, 0.157590740144327)
            ],
            (0.0, 1.0): [
                (0.0, 0.61573709),
                (0.87006477, 0.12993523)
            ],
            (0.98341909, 0.0): [
                (0.0, 0.0),
                (0.731724423625055, 0.157590737804387),
                (1.0, 0.0)
            ],
            (1.0, 0.0): [
                (0.88556116, 0.11443884),
                (0.98341909, 0.0)
            ],
            (0.731724419887841, 0.157590740144327): [
                (0.0, 0.36284274), (0.0, 0.61573709),
                (0.731724423625055, 0.157590737804387),
                (0.731724439651225, 0.157590734600608)
            ],
            (0.731724423625055, 0.157590737804387): [
                (0.0, 0.30386916),
                (0.731724419887841, 0.157590740144327),
                (0.731724439651225, 0.157590734600608),
                (0.98341909, 0.0)
            ],
            (0.731724439651225, 0.157590734600608): [
                (0.731724419887841, 0.157590740144327),
                (0.731724423625055, 0.157590737804387),
                (0.87006477, 0.12993523),
                (0.88556116, 0.11443884)
            ],
            (0.87006477, 0.12993523): [
                (0.0, 1.0),
                (0.731724439651225, 0.157590734600608),
                (0.88556116, 0.11443884)
            ],
            (0.88556116, 0.11443884): [
                (0.731724439651225, 0.157590734600608),
                (0.87006477, 0.12993523),
                (1.0, 0.0)
            ]
        },
        'search_edges': [
            (
                (0.0, 0.0),
                (0.0, 0.30386916)
            ),
            (
                (0.0, 0.0),
                (0.98341909, 0.0)
            ),
            (
                (0.0, 0.30386916),
                (0.0, 0.36284274)
            ),
            (
                (0.0, 0.30386916),
                (0.731724423625055, 0.157590737804387)
            ),
            (
                (0.0, 0.30386916),
                (0.731724423625055, 0.157590737804387)
            ),
            (
                (0.0, 0.36284274),
                (0.0, 0.61573709)
            ),
            (
                (0.0, 0.36284274),
                (0.731724419887841, 0.157590740144327)),
            (
                (0.0, 0.36284274),
                (0.731724419887841, 0.157590740144327)
            ),
            (
                (0.0, 0.61573709),
                (0.0, 1.0)
            ),
            (
                (0.0, 0.61573709),
                (0.731724419887841, 0.157590740144327)
            ),
            (
                (0.0, 0.61573709),
                (0.731724419887841, 0.157590740144327)
            ),
            (
                (0.0, 1.0),
                (0.87006477, 0.12993523)
            ),
            (
                (0.731724419887841, 0.157590740144327),
                (0.731724423625055, 0.157590737804387)
            ),
            (
                (0.731724419887841, 0.157590740144327),
                (0.731724423625055, 0.157590737804387)
            ),
            (
                (0.731724419887841, 0.157590740144327),
                (0.731724439651225, 0.157590734600608)
            ),
            (
                (0.731724419887841, 0.157590740144327),
                (0.731724439651225, 0.157590734600608)
            ),
            (
                (0.731724423625055, 0.157590737804387),
                (0.731724439651225, 0.157590734600608)),
            (
                (0.731724423625055, 0.157590737804387),
                (0.731724439651225, 0.157590734600608)),
            (
                (0.731724423625055, 0.157590737804387),
                (0.98341909, 0.0)),
            (
                (0.731724423625055, 0.157590737804387),
                (0.98341909, 0.0)),
            (
                (0.731724439651225, 0.157590734600608),
                (0.87006477, 0.12993523)),
            (
                (0.731724439651225, 0.157590734600608),
                (0.87006477, 0.12993523)),
            (
                (0.731724439651225, 0.157590734600608),
                (0.88556116, 0.11443884)),
            (
                (0.731724439651225, 0.157590734600608),
                (0.88556116, 0.11443884)),
            (
                (0.87006477, 0.12993523),
                (0.88556116, 0.11443884)
            ),
            (
                (0.88556116, 0.11443884),
                (1.0, 0.0)
            ),
            (
                (0.98341909, 0.0),
                (1.0, 0.0)
            )
        ],
        'search_nodes': [
            (0.0, 0.0),
            (0.0, 0.30386916),
            (0.0, 0.30386916),
            (0.0, 0.36284274),
            (0.0, 0.36284274),
            (0.0, 0.61573709),
            (0.0, 0.61573709),
            (0.0, 1.0),
            (0.731724419887841, 0.157590740144327),
            (0.731724419887841, 0.157590740144327),
            (0.731724419887841, 0.157590740144327),
            (0.731724419887841, 0.157590740144327),
            (0.731724423625055, 0.157590737804387),
            (0.731724423625055, 0.157590737804387),
            (0.731724423625055, 0.157590737804387),
            (0.731724423625055, 0.157590737804387),
            (0.731724439651225, 0.157590734600608),
            (0.731724439651225, 0.157590734600608),
            (0.731724439651225, 0.157590734600608),
            (0.731724439651225, 0.157590734600608),
            (0.87006477, 0.12993523),
            (0.87006477, 0.12993523),
            (0.88556116, 0.11443884),
            (0.88556116, 0.11443884),
            (0.98341909, 0.0),
            (0.98341909, 0.0),
            (1.0, 0.0)
        ],
        'polygons': [
            [
                (0.0, 0.0),
                (0.0, 0.30386916),
                (0.731724423625055, 0.157590737804387),
                (0.98341909, 0.0),
                (0.0, 0.0)
            ],
            [
                (0.0, 0.30386916),
                (0.0, 0.36284274),
                (0.731724419887841, 0.157590740144327),
                (0.731724423625055, 0.157590737804387),
                (0.0, 0.30386916)
            ],
            [
                (0.0, 0.36284274),
                (0.0, 0.61573709),
                (0.731724419887841, 0.157590740144327),
                (0.0, 0.36284274)
            ],
            [
                (0.0, 0.61573709),
                (0.0, 1.0),
                (0.87006477, 0.12993523),
                (0.731724439651225, 0.157590734600608),
                (0.731724419887841, 0.157590740144327),
                (0.0, 0.61573709)
            ],
            [
                (0.731724419887841, 0.157590740144327),
                (0.731724423625055, 0.157590737804387),
                (0.731724439651225, 0.157590734600608),
                (0.731724419887841, 0.157590740144327)
            ],
            [
                (0.731724423625055, 0.157590737804387),
                (0.731724439651225, 0.157590734600608),
                (0.88556116, 0.11443884),
                (1.0, 0.0),
                (0.98341909, 0.0),
                (0.731724423625055, 0.157590737804387)
            ],
            [
                (0.731724439651225, 0.157590734600608),
                (0.87006477, 0.12993523),
                (0.88556116, 0.11443884),
                (0.731724439651225, 0.157590734600608)
            ]
        ],
        'predictor_areas': {
            'PolyPhen-2': 0.11410069066319463,
            'SIFT': 0.1895963750902457,
            'CADD': 0.19630293424655962
        },
        'predictor_relative_areas': {
            'PolyPhen-2': 0.22820138132638926,
            'SIFT': 0.3791927501804914,
            'CADD': 0.39260586849311924
        },
        'output': '''
CLINICAL SPACE PARTITION
------------------------

Parameters considered: sensitivity, specificity and coverage

Methods compared: PolyPhen-2, SIFT, CADD

Best combination of methods (rho=0.5): CADD, SIFT, PolyPhen-2

List of clinical space fraction for each predictor:

Predictor 	Absolute value	Relative value
--------- 	--------------	--------------
CADD      	0.196		0.393
SIFT      	0.19		0.379
PolyPhen-2	0.114		0.228
'''

    }


def test_parse_config(base_case):
    rho, predictors = csp_co.parse_config(base_case['filename'], mode='co')
    assert rho == base_case['rho']
    assert predictors == base_case['predictors_demo']


def test_main(base_case, capsys):
    csp_co.main(base_case['rho'], base_case['predictors'])
    captured = capsys.readouterr()
    assert captured.out == base_case['output']


def test_get_predictors_intersection(base_case):
    nodes, lines = find_predictor_intersections.get_predictors_intersection(base_case['rho'], base_case['predictors'],
                                                                            precision=8)
    assert nodes == base_case['nodes']
    assert lines == base_case['lines']


def test_get_predictors_graph(base_case):
    interactions, search_edges, search_nodes = build_intersection_graph.get_predictors_graph(base_case['lines'])
    assert interactions == base_case['interactions']
    assert search_edges == base_case['search_edges']
    assert search_nodes == base_case['search_nodes']


def test_get_polygons(base_case):
    polygons = search_graph_polygons.get_polygons(base_case['search_edges'], base_case['search_nodes'],
                                                  base_case['nodes'], base_case['interactions'])
    assert polygons == base_case['polygons']


def test_get_polygons_data(base_case):
    predictor_areas, predictor_relative_areas = obtain_polygon_data.get_polygons_data(
        base_case['rho'], base_case['predictors'], base_case['polygons'])
    assert predictor_areas == base_case['predictor_areas']
    assert predictor_relative_areas == base_case['predictor_relative_areas']