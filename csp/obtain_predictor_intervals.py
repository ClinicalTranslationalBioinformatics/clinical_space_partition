"""
Get the intersection of predictors, the best predictor in each interval and merge intervals with same best predictor
"""

import numpy as np
from collections import defaultdict
from itertools import combinations


def get_intersection_point(x_1, y_1, x_2, y_2):
    """
    Get the intersection point of two predictors
    """
    m_1 = np.array([[x_1, -1], [x_2, -1]])
    n_1 = np.array([-y_1, -y_2])
    m_2 = np.array([[x_2, -1], [x_1, -1]])
    n_2 = np.array([-y_2, -y_1])
    try:
        x_a = np.linalg.solve(m_1, n_1).tolist()[0]
        x_b = np.linalg.solve(m_2, n_2).tolist()[0]
    except np.linalg.LinAlgError as e:
        if str(e) == 'Singular matrix':
            return None
    x = (x_a + x_b) / 2
    if 0 <= x <= 1:
        return x
    else:
        return None


def get_predictors_intersections(rho, predictors):
    """
    Get the intersection point between each pair of predictors
    """
    x_points = set()
    for predictor1, predictor2 in combinations(predictors, 2):
        sens_1, spec_1 = predictors[predictor1]
        sens_2, spec_2 = predictors[predictor2]
        x_1 = 1 - spec_1 + rho * (sens_1 + spec_1 - 2)
        y_1 = rho * (1 - sens_1)
        x_2 = 1 - spec_2 + rho * (sens_2 + spec_2 - 2)
        y_2 = rho * (1 - sens_2)
        x = get_intersection_point(x_1, y_1, x_2, y_2)
        if x:
            x_points.add(x)
    return sorted(x_points)


def get_middle_points(x_points):
    """
    Get the middle point of an interval
    """
    middle_points = {}
    if len(x_points) == 0:
        middle_points[(0, 1)] = 0.5
    elif len(x_points) == 1:
        value = x_points[0]
        middle_points[(0, value)] = value / 2
        middle_points[(value, 1)] = (value + 1) / 2
    else:
        for index, value in enumerate(x_points):
            if index + 1 == len(x_points):
                middle_points[(value, 1)] = (value + 1) / 2
            else:
                next_value = x_points[index + 1]
                middle_points[(value, next_value)] = (value + next_value) / 2
                if index == 0:
                    middle_points[(0, value)] = value / 2
    return middle_points


def get_interval_best_predictor(rho, predictors, x_points):
    """
    Get the best predictor of each interval
    """
    middle_points = get_middle_points(x_points)
    interval_best_predictors = {}
    for interval, middle_point in middle_points.items():
        cost_predictors = {predictor: middle_point * ((1 - spec) + rho * (sens + spec - 2)) + rho * (1 - sens)
                           for predictor, (sens, spec) in predictors.items()}
        interval_best_predictors[interval] = min(cost_predictors, key=cost_predictors.get)
    return interval_best_predictors


def merge_intervals(interval_best_predictors):
    """
    Merge intervals with the same best predictor
    """
    predictor2intervals = defaultdict(set)
    for interval, best_predictor in interval_best_predictors.items():
        predictor2intervals[best_predictor].update(interval)
    merged_intervals = {best_predictor: max(interval_points) - min(interval_points)
                        for best_predictor, interval_points in predictor2intervals.items()}
    return merged_intervals
