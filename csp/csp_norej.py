"""
Cost space partition without coverage
"""

from csp_rej import parse_args, parse_config, print_float
from obtain_predictor_intervals import get_predictors_intersections, get_interval_best_predictor, merge_intervals


def print_output(rho, predictors, merged_intervals):
    merged_intervals.update({predictor: 0 for predictor in predictors if predictor not in merged_intervals})
    sorted_intervals = sorted(merged_intervals.items(), key=lambda x: (-x[1], x[0]))
    predictors_area_round3 = [predictor for predictor, interval in sorted_intervals if round(interval, 3) > 0]
    spaces_predictors = len(max(list(predictors) + ['Predictor'], key=lambda p: len(p)))
    print('\nCLINICAL SPACE PARTITION')
    print('------------------------\n')
    print('Parameters considered: sensitivity and specificity\n')
    print('Methods compared: {}\n'.format(', '.join(predictors)))
    print('Best combination of methods (rho={}): {}\n'.format(rho, ', '.join(predictors_area_round3)))
    print('List of clinical space fraction for each predictor:\n')
    print('{: <{spaces}}\tRelative value'.format('Predictor', spaces=spaces_predictors))
    print('{: <{spaces}}\t--------------'.format('---------', spaces=spaces_predictors))
    for best_predictor, merged_interval in sorted_intervals:
        print('{: <{spaces}}\t{}'.format(best_predictor, print_float(merged_interval), spaces=spaces_predictors))


def main(rho, predictors):
    """
    Get the cost space partition of given predictors without coverage
    """
    # Get the intersection points of predictors
    x_points = get_predictors_intersections(rho, predictors)

    # Get the best predictor in each interval
    interval_best_predictors = get_interval_best_predictor(rho, predictors, x_points)

    # Merge intervals with the same best predictor
    merged_intervals = merge_intervals(interval_best_predictors)

    # Output
    print_output(rho, predictors, merged_intervals)


if __name__ == '__main__':
    # Parse predictors and rho
    user_rho, user_predictors = parse_config(parse_args(), mode='norej')

    # Execute CSP without coverage
    main(user_rho, user_predictors)
