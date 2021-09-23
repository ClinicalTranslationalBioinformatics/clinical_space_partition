"""
Cost space partition with coverage
"""

import sys
import argparse
import configparser
from find_predictor_intersections import get_predictors_intersection
from build_intersection_graph import get_predictors_graph
from search_graph_polygons import get_polygons
from obtain_polygon_data import get_polygons_data
import decimal


def parse_args():
    """
    Parse command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=argparse.FileType('r'), help='select the config file')
    args = parser.parse_args()
    filename = args.file.name
    return filename


def parse_config(filename, mode):
    """
    Parse the rho and predictor's sensitivity, specificity and coverage from the config file
    """
    config = configparser.ConfigParser()
    config.optionxform = str
    try:
        config.read(filename)
        rho = float(config.get('rho', 'rho'))
        predictors = dict(config.items('predictors'))
    except Exception as e:
        sys.exit(e)

    if not 0.00001 <= rho <= 1:
        sys.exit('The rho value ' + str(rho) + ' should be between 0.00001 - 1')

    if len(predictors) == 0:
        sys.exit('The predictor(s) are missing')
    if mode == 'co':
        info = ['sensitivity', 'specificity', 'coverage']
    elif mode == 'noco':
        info = ['sensitivity', 'specificity']
    else:
        raise Exception('ERROR: config mode unknown')
    for predictor, values in predictors.items():
        values = values.replace(' ', '').split(',')
        if len(values) != len(info):
            message_values = ', '.join(info[:-1]) + ' and ' + info[-1]
            sys.exit(predictor + ' information should contain ' + message_values + ' but it has ' + str(len(values)) +
                     ' elements')
        for i, value in enumerate(values):
            try:
                float(value)
            except ValueError:
                sys.exit(predictor + ' ' + info[i] + ' is ' + str(value) + ' but should be a number')
        values = [round(float(value), 3) for value in values]
        for i, value in enumerate(values):
            if not 0 <= value <= 1:
                sys.exit(predictor + ' ' + info[i] + ' is ' + str(value) + ' but should be between 0 - 1')
        predictors[predictor] = values

    return rho, predictors


def predictors_2_polygons(rho, predictors, precision):
    """
    Get the polygons from the intersection of predictors
    """
    # Get the lines and nodes of the intersection of predictor's planes and lines
    nodes, lines = get_predictors_intersection(rho, predictors, precision)

    # Get graph elements to search the polygons
    interactions, search_edges, search_nodes = get_predictors_graph(lines)

    # Search the polygons
    polygons = get_polygons(search_edges, search_nodes, nodes, interactions)

    return polygons


def print_float(num):
    """
    Get float with 3 decimals and normalized
    """
    return decimal.Decimal(str(round(num, 3))).normalize()


def print_output(rho, predictors, predictor_areas, predictor_relative_areas):
    sorted_predictor_areas = sorted(predictor_relative_areas.items(), key=lambda x: (-x[1], x[0]))
    predictors_area_round3 = [predictor for predictor, area in sorted_predictor_areas if round(area, 3) > 0]
    spaces_predictors = len(max(list(predictors) + ['Predictor'], key=lambda p: len(p)))
    print('\nCLINICAL SPACE PARTITION')
    print('------------------------\n')
    print('Parameters considered: sensitivity, specificity and coverage\n')
    print('Methods compared: {}\n'.format(', '.join(predictors)))
    print('Best combination of methods (rho={}): {}\n'.format(rho, ', '.join(predictors_area_round3)))
    print('List of clinical space fraction for each predictor:\n')
    print('{: <{spaces}}\tAbsolute value\tRelative value'.format('Predictor', spaces=spaces_predictors))
    print('{: <{spaces}}\t--------------\t--------------'.format('---------', spaces=spaces_predictors))
    for predictor, _ in sorted_predictor_areas:
        print('{: <{spaces}}\t{}\t\t{}'.format(predictor, print_float(predictor_areas.get(predictor, 0)),
                                               print_float(predictor_relative_areas.get(predictor, 0)),
                                               spaces=spaces_predictors))


def main(rho, predictors):
    """
    Get the cost space partition of given predictors with coverage
    """
    try:
        polygons = predictors_2_polygons(rho, predictors, precision=8)
    except IndexError:
        polygons = predictors_2_polygons(rho, predictors, precision=10)

    # Get best predictors, areas and relative areas
    predictor_areas, predictor_relative_areas = get_polygons_data(rho, predictors, polygons)

    # Output
    print_output(rho, predictors, predictor_areas, predictor_relative_areas)


if __name__ == '__main__':
    # Parse predictors and rho
    user_rho, user_predictors = parse_config(parse_args(), mode='co')

    # Execute CSP coverage
    main(user_rho, user_predictors)
