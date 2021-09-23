"""
Calculate the middle points, areas and best predictors of the polygons and the relative areas of the best predictors
"""

from shapely.geometry.polygon import Polygon


def get_predictor_cost(x, y, rho, sens, spec, cov):
    """
    Calculate the predictor's cost on a point (x, y) based on the rho and its sensitivity, specificity and coverage
    """
    return x * ((rho * cov * (1 - sens)) + cov - 1) + y * (((1 - rho) * cov * (1 - spec)) + cov - 1) + 1 - cov


def get_polygon_best_predictor(rho, predictors, polygons):
    """
    Calculate the predictor with the best cost in a polygon
    """
    polygon_best_predictor = []
    for polygon_id, polygon in enumerate(polygons):
        centroid = Polygon(polygon).centroid
        cost_predictors = {predictor: get_predictor_cost(centroid.x, centroid.y, rho, *predictors[predictor])
                           for predictor in predictors}
        best_predictor = min(cost_predictors, key=cost_predictors.get)
        polygon_best_predictor.append(best_predictor)
    return polygon_best_predictor


def get_best_predictor_area(polygons, polygon_best_predictor):
    """
    Calculate the total area of the best predictors in the triangle
    """
    best_predictor_areas = dict.fromkeys(polygon_best_predictor, 0.0)
    for polygon_id, polygon in enumerate(polygons):
        best_predictor_areas[polygon_best_predictor[polygon_id]] += Polygon(polygon).area
    return best_predictor_areas


def get_predictor_area(best_predictor_areas, predictors):
    """
    Calculate the area and relative areas of predictors
    """
    predictor_areas = {predictor: best_predictor_areas.get(predictor, 0.0) for predictor in predictors}
    predictor_relative_areas = {predictor: predictor_areas[predictor] / 0.5 for predictor in predictors}
    return predictor_areas, predictor_relative_areas


def get_polygons_data(rho, predictors, polygons):
    """
    Calculate the best predictor in a polygon and best predictors' areas and relative areas
    """
    polygon_best_predictor = get_polygon_best_predictor(rho, predictors, polygons)
    best_predictor_areas = get_best_predictor_area(polygons, polygon_best_predictor)
    predictor_areas, predictor_relative_areas = get_predictor_area(best_predictor_areas, predictors)
    return predictor_areas, predictor_relative_areas
