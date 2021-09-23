from csp import csp_co, find_predictor_intersections, search_graph_polygons


def test_linear_equation_divided_zero_values():
    slope, intercept = find_predictor_intersections.get_linear_equation_parameters(0.5, 0.965, 0.76, 1, 0.899,
                                                                                   0.760, 1)
    assert slope is None
    assert intercept is None


def test_line_merges_x_axis():
    rho = 0.5
    predictors = {
        'predictor1': [0.0, 0.0, 0.0],
        'predictor2': [0.0, 0.0, 0.003],
        'predictor3': [0.0, 0.003, 0.003]
    }
    expected_lines = {
        'hypotenuse': [
            (0.0, 1.0),
            (1.0, 0.0)
        ],
        'x_axis': [
            (0.0, 0.0),
            (0.666666666666667, 0.0),
            (0.66666667, 0.0),
            (1.0, 0.0)
        ],
        'y_axis': [
            (0.0, 0.0),
            (0.0, 0.66666667),
            (0.0, 0.667334),
            (0.0, 1.0)
        ],
        (-1.001001, 0.667334): [
            (0.0, 0.667334),
            (0.666663336663285, 3.333336715e-06),
            (0.666666666666667, 0.0)],
        (-1.0, 0.66666667): [
            (0.0, 0.66666667),
            (0.666663336663285, 3.333336715e-06),
            (0.66666667, 0.0)]
    }
    nodes, lines = find_predictor_intersections.get_predictors_intersection(rho, predictors, precision=8)
    assert lines == expected_lines


def test_line_merges_hypotenuse():
    rho = 0.5
    predictors = {
        'predictor1': [0.0, 0.0, 0.0],
        'predictor2': [1.0, 1.0, 1.0],
        'predictor3': [0.0, 0.0, 0.5]
    }
    expected_lines = {
        'y_axis': [
            (0.0, 0.0),
            (0.0, 0.66666667),
            (0.0, 1.0)
        ],
        'x_axis': [
            (0.0, 0.0),
            (0.66666667, 0.0),
            (1.0, 0.0)
        ],
        'hypotenuse': [
            (0.0, 1.0),
            (1.0, 0.0)
        ],
        (-1.0, 0.66666667): [
            (0.0, 0.66666667),
            (0.66666667, 0.0)
        ]
    }
    nodes, lines = find_predictor_intersections.get_predictors_intersection(rho, predictors, precision=8)
    assert lines == expected_lines


def test_perfect_predictor(capsys):
    rho = 0.5
    predictors = {
        'predictor1': [1.0, 1.0, 1.0],
        'predictor2': [0.83, 0.92, 1.0],
        'predictor3': [0.95, 0.95, 0.9]
    }
    output = '''
CLINICAL SPACE PARTITION
------------------------

Parameters considered: sensitivity, specificity and coverage

Methods compared: predictor1, predictor2, predictor3

Best combination of methods (rho=0.5): predictor1

List of clinical space fraction for each predictor:

Predictor 	Absolute value	Relative value
--------- 	--------------	--------------
predictor1	0.5		1
predictor2	0		0
predictor3	0		0
'''
    csp_co.main(rho, predictors)
    captured = capsys.readouterr()
    assert captured.out == output


def test_worst_predictors(capsys):
    rho = 0.5
    predictors = {
        'predictor1': [0, 0, 0],
        'predictor2': [0, 0, 0],
        'predictor3': [0, 0, 0]
    }
    output = '''
CLINICAL SPACE PARTITION
------------------------

Parameters considered: sensitivity, specificity and coverage

Methods compared: predictor1, predictor2, predictor3

Best combination of methods (rho=0.5): predictor1

List of clinical space fraction for each predictor:

Predictor 	Absolute value	Relative value
--------- 	--------------	--------------
predictor1	0.5		1
predictor2	0		0
predictor3	0		0
'''
    csp_co.main(rho, predictors)
    captured = capsys.readouterr()
    assert captured.out == output


def test_search_polygons_c5(capsys):
    rho = 0.00001
    predictors = {
        'CADD': [0.995, 0.254, 1.0],
        'MetaLR': [0.887, 0.853, 0.996],
        'MutationTaster': [0.979, 0.603, 0.996]
    }
    output = '''
CLINICAL SPACE PARTITION
------------------------

Parameters considered: sensitivity, specificity and coverage

Methods compared: CADD, MetaLR, MutationTaster

Best combination of methods (rho=1e-05): MetaLR, CADD

List of clinical space fraction for each predictor:

Predictor     	Absolute value	Relative value
---------     	--------------	--------------
MetaLR        	0.497		0.993
CADD          	0.003		0.007
MutationTaster	0		0
'''
    csp_co.main(rho, predictors)
    captured = capsys.readouterr()
    assert captured.out == output


def test_search_polygons_c6():
    nodes = {
        (0.0, 0.0): {'x_axis', 'y_axis'},
        (0.0, 1.0): {'y_axis', 'hypotenuse'},
        (1.0, 0.0): {'x_axis', 'hypotenuse'},
        (0.0, 0.4): {'y_axis', 'y=0.4'},
        (0.0, 0.6): {'y_axis', 'y=-x+0.6', 'y=0.6'},
        (0.2, 0.0): {'x_axis', 'x=0.2', 'y=x-0.2'},
        (0.6, 0.0): {'x_axis', 'x=0.6', 'y=-x+0.6'},
        (0.2, 0.8): {'hypotenuse', 'x=0.2'},
        (0.6, 0.4): {'hypotenuse', 'x=0.6', 'y=x-0.2', 'y=0.4'},
        (0.2, 0.4): {'x=0.2', 'y=-x+0.6', 'y=0.4'},
        (0.4, 0.2): {'y=-x+0.6', 'y=x-0.2'},
        (0.4, 0.6): {'hypotenuse', 'y=0.6'},
        (0.2, 0.6): {'x=0.2', 'y=0.6'}
    }

    polygon_edges = [
        ((0.0, 0.0), (0.0, 0.4)), ((0.0, 0.0), (0.2, 0.0)), ((0.0, 0.4), (0.0, 0.6)), ((0.0, 0.4), (0.2, 0.4)),
        ((0.0, 0.4), (0.2, 0.4)), ((0.0, 0.6), (0.0, 1.0)), ((0.0, 0.6), (0.2, 0.4)), ((0.0, 0.6), (0.2, 0.4)),
        ((0.0, 0.6), (0.2, 0.6)), ((0.0, 0.6), (0.2, 0.6)), ((0.0, 1.0), (0.2, 0.8)), ((0.2, 0.0), (0.2, 0.4)),
        ((0.2, 0.0), (0.2, 0.4)), ((0.2, 0.0), (0.4, 0.2)), ((0.2, 0.0), (0.4, 0.2)), ((0.2, 0.0), (0.6, 0.0)),
        ((0.2, 0.4), (0.2, 0.6)), ((0.2, 0.4), (0.2, 0.6)), ((0.2, 0.4), (0.4, 0.2)), ((0.2, 0.4), (0.4, 0.2)),
        ((0.2, 0.4), (0.6, 0.4)), ((0.2, 0.4), (0.6, 0.4)), ((0.2, 0.6), (0.2, 0.8)), ((0.2, 0.6), (0.2, 0.8)),
        ((0.2, 0.6), (0.4, 0.6)), ((0.2, 0.6), (0.4, 0.6)), ((0.2, 0.8), (0.4, 0.6)), ((0.4, 0.2), (0.6, 0.0)),
        ((0.4, 0.2), (0.6, 0.0)), ((0.4, 0.2), (0.6, 0.4)), ((0.4, 0.2), (0.6, 0.4)), ((0.4, 0.6), (0.6, 0.4)),
        ((0.6, 0.0), (0.6, 0.4)), ((0.6, 0.0), (0.6, 0.4)), ((0.6, 0.0), (1.0, 0.0)), ((0.6, 0.4), (1.0, 0.0))
    ]

    interactions = {
        (0.0, 0.0): [(0.2, 0.0), (0.0, 0.4)],
        (0.0, 1.0): [(0.0, 0.6), (0.2, 0.8)],
        (1.0, 0.0): [(0.6, 0.0), (0.6, 0.4)],
        (0.0, 0.4): [(0.0, 0.0), (0.2, 0.4), (0.0, 0.6)],
        (0.0, 0.6): [(0.0, 0.4), (0.2, 0.4), (0.2, 0.6), (0.0, 1.0)],
        (0.2, 0.0): [(0.0, 0.0), (0.2, 0.4), (0.4, 0.2), (0.6, 0.0)],
        (0.6, 0.0): [(0.2, 0.0), (0.4, 0.2), (0.6, 0.4), (1.0, 0.0)],
        (0.2, 0.8): [(0.0, 1.0), (0.2, 0.6), (0.4, 0.6)],
        (0.6, 0.4): [(0.6, 0.0), (0.4, 0.6), (0.2, 0.4), (0.4, 0.2), (1.0, 0.0)],
        (0.2, 0.4): [(0.6, 0.4), (0.2, 0.6), (0.0, 0.6), (0.0, 0.4), (0.2, 0.0), (0.4, 0.2)],
        (0.4, 0.2): [(0.2, 0.4), (0.2, 0.0), (0.6, 0.0), (0.6, 0.4)],
        (0.4, 0.6): [(0.2, 0.8), (0.2, 0.6), (0.6, 0.4)],
        (0.2, 0.6): [(0.2, 0.8), (0.4, 0.6), (0.2, 0.4), (0.0, 0.6)]
    }

    polygon_nodes = [
        (0.2, 0.4), (0.2, 0.4), (0.2, 0.4), (0.2, 0.4), (0.2, 0.4), (0.2, 0.4), (0.0, 0.0), (0.0, 1.0),
        (1.0, 0.0), (0.0, 0.4), (0.0, 0.4), (0.0, 0.6), (0.0, 0.6), (0.0, 0.6), (0.2, 0.0), (0.2, 0.0),
        (0.2, 0.0), (0.6, 0.0), (0.6, 0.0), (0.6, 0.0), (0.2, 0.8), (0.2, 0.8), (0.6, 0.4), (0.6, 0.4),
        (0.6, 0.4), (0.6, 0.4), (0.4, 0.2), (0.4, 0.2), (0.4, 0.2), (0.4, 0.2), (0.2, 0.6), (0.2, 0.6),
        (0.2, 0.6), (0.2, 0.6), (0.4, 0.6), (0.4, 0.6)
    ]

    expected_polygons = [
        [(0.2, 0.4), (0.6, 0.4), (0.4, 0.2), (0.2, 0.4)],
        [(0.2, 0.4), (0.2, 0.6), (0.0, 0.6), (0.2, 0.4)],
        [(0.2, 0.4), (0.0, 0.6), (0.0, 0.4), (0.2, 0.4)],
        [(0.2, 0.4), (0.2, 0.0), (0.4, 0.2), (0.2, 0.4)],
        [(0.2, 0.4), (0.6, 0.4), (0.4, 0.6), (0.2, 0.6), (0.2, 0.4)],
        [(0.2, 0.4), (0.0, 0.4), (0.0, 0.0), (0.2, 0.0), (0.2, 0.4)],
        [(0.0, 1.0), (0.0, 0.6), (0.2, 0.6), (0.2, 0.8), (0.0, 1.0)],
        [(1.0, 0.0), (0.6, 0.0), (0.6, 0.4), (1.0, 0.0)],
        [(0.2, 0.0), (0.4, 0.2), (0.6, 0.0), (0.2, 0.0)],
        [(0.6, 0.0), (0.4, 0.2), (0.6, 0.4), (0.6, 0.0)],
        [(0.2, 0.8), (0.2, 0.6), (0.4, 0.6), (0.2, 0.8)]
    ]

    polygons = search_graph_polygons.get_polygons(polygon_edges, polygon_nodes, nodes, interactions)
    assert polygons == expected_polygons


def test_sort_line_y_axis(capsys):
    rho = 0.00001
    predictors = {
        'CADD': [0.995, 0.254, 1.0],
        'fathmm': [0.835, 0.658, 0.902],
        'MetaLR': [0.887, 0.853, 0.996],
        'MutPred': [0.95, 0.706, 0.281],
        'PolyPhen-2': [0.926, 0.638, 0.909],
        'SIFT': [0.924, 0.682, 0.866]
    }
    output = '''
CLINICAL SPACE PARTITION
------------------------

Parameters considered: sensitivity, specificity and coverage

Methods compared: CADD, fathmm, MetaLR, MutPred, PolyPhen-2, SIFT

Best combination of methods (rho=1e-05): MetaLR, MutPred, CADD

List of clinical space fraction for each predictor:

Predictor 	Absolute value	Relative value
--------- 	--------------	--------------
MetaLR    	0.456		0.911
MutPred   	0.041		0.082
CADD      	0.003		0.007
SIFT      	0		0
PolyPhen-2	0		0
fathmm    	0		0
'''
    csp_co.main(rho, predictors)
    captured = capsys.readouterr()
    assert captured.out == output


def test_black_hole(capsys):
    rho = 0.00001
    predictors = {
        'CADD': [0.995, 0.254, 1.0],
        'MetaLR': [0.887, 0.853, 0.996],
        'MetaSVM': [0.896, 0.879, 0.996],
        'MutPred': [0.95, 0.706, 0.281],
        'MutationTaster': [0.979, 0.603, 0.996],
        'PolyPhen-2': [0.926, 0.638, 0.909],
        'SIFT': [0.924, 0.682, 0.866]
    }
    output = '''
CLINICAL SPACE PARTITION
------------------------

Parameters considered: sensitivity, specificity and coverage

Methods compared: CADD, MetaLR, MetaSVM, MutPred, MutationTaster, PolyPhen-2, SIFT

Best combination of methods (rho=1e-05): MetaSVM, MutPred, CADD

List of clinical space fraction for each predictor:

Predictor     	Absolute value	Relative value
---------     	--------------	--------------
MetaSVM       	0.472		0.943
MutPred       	0.025		0.05
CADD          	0.003		0.006
MutationTaster	0		0
MetaLR        	0		0
PolyPhen-2    	0		0
SIFT          	0		0
'''
    csp_co.main(rho, predictors)
    captured = capsys.readouterr()
    assert captured.out == output


def test_dynamic_precision(capsys):
    rho = 0.00001
    predictors = {
        'predictor1': [0.163, 0.884, 0.724],
        'predictor2': [0.242, 0.176, 0.828],
        'predictor3': [0.321, 0.946, 0.948],
        'predictor4': [0.093, 0.884, 0.67]
    }
    output = '''
CLINICAL SPACE PARTITION
------------------------

Parameters considered: sensitivity, specificity and coverage

Methods compared: predictor1, predictor2, predictor3, predictor4

Best combination of methods (rho=1e-05): predictor3

List of clinical space fraction for each predictor:

Predictor 	Absolute value	Relative value
--------- 	--------------	--------------
predictor3	0.5		1
predictor1	0		0
predictor4	0		0
predictor2	0		0
'''
    csp_co.main(rho, predictors)
    captured = capsys.readouterr()
    assert captured.out == output


def test_merge_nodes(capsys):
    rho = 0.00001
    predictors = {
        'predictor1': [0.111, 0.655, 0.098],
        'predictor2': [0.54, 0.92, 0.401],
        'predictor3': [0.328, 0.389, 0.883],
        'predictor4': [0.559, 0.751, 0.54],
        'predictor5': [0.38, 0.345, 0.382],
        'predictor6': [0.02, 0.304, 0.243]
    }
    output = '''
CLINICAL SPACE PARTITION
------------------------

Parameters considered: sensitivity, specificity and coverage

Methods compared: predictor1, predictor2, predictor3, predictor4, predictor5, predictor6

Best combination of methods (rho=1e-05): predictor3, predictor2, predictor4

List of clinical space fraction for each predictor:

Predictor 	Absolute value	Relative value
--------- 	--------------	--------------
predictor3	0.229		0.459
predictor2	0.212		0.424
predictor4	0.059		0.117
predictor1	0		0
predictor5	0		0
predictor6	0		0
'''
    csp_co.main(rho, predictors)
    captured = capsys.readouterr()
    assert captured.out == output


def test_merge_recursive_nodes(capsys):
    rho = 0.00001
    predictors = {
        'predictor1': [0.503, 0.846, 0.632],
        'predictor2': [0.584, 0.973, 0.522],
        'predictor3': [0.424, 0.101, 0.9],
        'predictor4': [0.504, 0.846, 0.638],
        'predictor5': [0.364, 0.843, 0.0],
        'predictor6': [0.813, 0.347, 0.982]
    }
    output = '''
CLINICAL SPACE PARTITION
------------------------

Parameters considered: sensitivity, specificity and coverage

Methods compared: predictor1, predictor2, predictor3, predictor4, predictor5, predictor6

Best combination of methods (rho=1e-05): predictor2, predictor6, predictor4, predictor5

List of clinical space fraction for each predictor:

Predictor 	Absolute value	Relative value
--------- 	--------------	--------------
predictor2	0.197		0.394
predictor6	0.194		0.388
predictor4	0.096		0.192
predictor5	0.013		0.026
predictor1	0		0
predictor3	0		0
'''
    csp_co.main(rho, predictors)
    captured = capsys.readouterr()
    assert captured.out == output


def test_sympy_lin_solve():
    rho = 0.00001
    predictors = {
        'predictor1': [0.942, 0.881, 0.996],
        'predictor2': [0.0, 0.0, 0.0],
        'predictor3': [0.003, 0.922, 0.218]
    }
    expected_polygons = [
        [
            (0.0, 0.0),
            (0.0, 0.8845744247),
            (0.999743205223076, 0.000228968026135),
            (0.999743205223091, 0.000228968026121),
            (0.9999900301, 0.0),
            (0.0, 0.0)
        ],
        [
            (0.0, 0.8845744247),
            (0.0, 0.8936559995),
            (0.999743205223076, 0.000228968026135),
            (0.0, 0.8845744247)
        ],
        [
            (0.0, 0.8936559995),
            (0.0, 0.927644456),
            (0.999743205223086, 0.000228968026126),
            (0.999743205223076, 0.000228968026135),
            (0.0, 0.8936559995)
        ],
        [
            (0.0, 0.927644456),
            (0.0, 1.0),
            (0.999984281, 1.5719e-05),
            (0.999743205223086, 0.000228968026126),
            (0.0, 0.927644456)
        ],
        [
            (0.999743205223076, 0.000228968026135),
            (0.999743205223086, 0.000228968026126),
            (0.999743205223091, 0.000228968026121),
            (0.999743205223076, 0.000228968026135)
        ],
        [
            (0.999743205223086, 0.000228968026126),
            (0.999743205223091, 0.000228968026121),
            (0.99999942, 0.0),
            (1.0, 0.0),
            (0.999984281, 1.5719e-05),
            (0.999743205223086, 0.000228968026126)
        ],
        [
            (0.999743205223091, 0.000228968026121),
            (0.9999900301, 0.0),
            (0.99999942, 0.0),
            (0.999743205223091, 0.000228968026121)
        ]
    ]
    polygons = csp_co.predictors_2_polygons(rho, predictors, 10)
    assert polygons == expected_polygons
