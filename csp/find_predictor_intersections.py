"""
Get the lines and nodes from the predictor's intersection
"""

from collections import defaultdict
from itertools import combinations
import math
import sympy

TRIANGLE_LINES = {'x_axis', 'hypotenuse', 'y_axis'}


def nodes2edge(n1, n2):
    """
    Get the edge of two nodes
    """
    return tuple(sorted([n1, n2], key=lambda x: (x[0], x[1])))


def get_linear_equation_parameters(rho, sens_1, spec_1, cov_1, sens_2, spec_2, cov_2):
    """
    Get the slope and intercept of the linear equation of the intersection of two planes (equation A3.3)
    ZeroDivisionError: planes are parallel
    """
    try:
        denominator = (cov_1 * (1 - rho) * (1 - spec_1) + cov_1 - 1) - (cov_2 * (1 - rho) * (1 - spec_2) + cov_2 - 1)
        slope = ((cov_2 * rho * (1 - sens_2) + cov_2 - 1) - (cov_1 * rho * (1 - sens_1) + cov_1 - 1)) / denominator
        intercept = ((1 - cov_2) - (1 - cov_1)) / denominator
        if str(intercept) == '-0.0':
            intercept = 0.0
        if str(slope) == '-0.0':
            slope = 0.0
        return slope, intercept
    except ZeroDivisionError:
        return None, None


def initialize_nodes():
    """
    Get the initial nodes and lines that form the triangle cost space
    """
    nodes = defaultdict(set, {(0.0, 0.0): {'x_axis', 'y_axis'},
                              (0.0, 1.0): {'y_axis', 'hypotenuse'},
                              (1.0, 0.0): {'hypotenuse', 'x_axis'}})
    return nodes


def get_intersection_point(slope1, intercept1, slope2, intercept2, precision):
    """
    Get the intersection point between two lines
    """
    slope = sympy.Matrix([[slope1, -1], [slope2, -1]])
    intercept = sympy.Matrix([-intercept1, -intercept2])
    resolution = sympy.linsolve((slope, intercept), [sympy.Symbol('x'), sympy.Symbol('y')])
    if resolution == sympy.EmptySet:
        return None
    x, y = list(resolution)[0]
    x, y = round(float(x), precision), round(float(y), precision)
    if not (0 <= x <= 1 and 0 <= y <= 1):
        return None
    if str(x) == '-0.0':
        x = 0.0
    if str(y) == '-0.0':
        y = 0.0
    point = (x, y)
    # Hypotenuse
    round_margin = 2 * 10 ** (-1 * precision)
    if slope2 == -1 and intercept2 == 1 and 1 - round_margin <= x + y <= 1 + round_margin:
        return point
    # Other lines
    elif slope2 != -1 and intercept2 != 1 and x + y <= 1:
        return point
    else:
        return None


def get_intersection_with_lines(slope1, intercept1, potential_lines, nodes, precision):
    """
    Get the intersection points between the linear equations defined by two planes
    and the set of linear equations defined by the rest of possible combinations of predictors
    """
    line1 = (slope1, intercept1)
    for line2 in potential_lines:
        if line1 == line2:
            continue
        slope2, intercept2 = line2
        point = get_intersection_point(slope1, intercept1, slope2, intercept2, precision)
        if point:
            nodes[point].update({line1, line2})
    return nodes


def get_intersection_point_xaxis(slope, intercept, precision):
    """
    Get the intersection point between a line and x axis
    """
    if slope != 0.0:
        x = round(-intercept / slope, precision)
        if str(x) == '-0.0':
            x = 0.0
        if 0.0 <= x <= 1.0:
            return (x, 0.0)
    return None


def get_intersection_point_yaxis(y):
    """
    Get the intersection point between a line and x axis
    """
    if str(y) == '-0.0':
        y = 0.0
    if 0.0 <= y <= 1.0:
        return (0.0, y)
    return None


def get_intersections(line, potential_lines, nodes, precision):
    """
    Get the intersection points between the lines defined by two planes
    and the lines defined by the cost area (x=0, y=0, y=-x+1)
    and the lines defined by the possible combinations of predictors
    """
    slope, intercept = line

    # Intersection with y axis when x=0
    point = get_intersection_point_yaxis(intercept)
    if point:
        nodes[point].update({line, 'y_axis'})

    # Intersection with x axis when y=0
    point = get_intersection_point_xaxis(slope, intercept, precision)
    if point:
        nodes[point].update({line, 'x_axis'})

    # Intersection with y=-x+1 line
    point = get_intersection_point(slope, intercept, -1, 1, precision)
    if point:
        nodes[point].update({line, 'hypotenuse'})

    # Intersection with other lines
    nodes = get_intersection_with_lines(slope, intercept, potential_lines, nodes, precision)
    return nodes


def get_nodes(rho, predictors, precision):
    """
    Get nodes
    """
    potential_lines = []
    for group in combinations(predictors, 2):
        slope, intercept = get_linear_equation_parameters(rho, *predictors[group[0]], *predictors[group[1]])
        if slope is not None and intercept is not None:
            line = (round(slope, precision), round(intercept, precision))
            if line in [(0.0, 0.0), (-1.0, 1.0)]:  # Discard lines equal to x_axis and hypotenuse
                continue
            elif line not in potential_lines:
                potential_lines.append(line)
    nodes = initialize_nodes()
    for line in potential_lines:
        if line in TRIANGLE_LINES:
            continue
        nodes = get_intersections(line, potential_lines, nodes, precision)
    return nodes


def unmerge_nodes(nodes):
    """
    Unmerge nodes with more than 2 lines with higher precision
    """
    lines2node = {tuple(sorted(map(str, lines))): node for node, lines in nodes.items()}
    for node, lines in dict(nodes).items():
        if len(lines) > 2:
            # No unmerge
            predictors_lines = lines - TRIANGLE_LINES
            if node == (0.0, 0.0) and all([line[1] == 0.0 for line in predictors_lines]):
                continue
            if 'y_axis' in lines and len(set([line[1] for line in predictors_lines])) == 1:
                continue
            # Unmerge
            del nodes[node]
            for lines_comb in combinations(sorted(lines, key=lambda line: str(line)), 2):
                # Remove other affected nodes
                lines_comb_id = tuple(sorted(map(str, lines_comb)))
                node_comb = lines2node.get(lines_comb_id)
                if lines_comb_id in lines2node and node_comb in nodes:
                    del nodes[node_comb]
                # Get new node
                if node in [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)] and set(lines_comb) - TRIANGLE_LINES == set():
                    point = node
                elif 'x_axis' in lines_comb:
                    other_line = [line for line in lines_comb if line not in TRIANGLE_LINES][0]
                    point = get_intersection_point_xaxis(other_line[0], other_line[1], precision=15)
                elif 'y_axis' in lines_comb:
                    other_line = [line for line in lines_comb if line not in TRIANGLE_LINES][0]
                    point = get_intersection_point_yaxis(other_line[1])
                elif 'hypotenuse' in lines_comb:
                    other_line = [line for line in lines_comb if line not in TRIANGLE_LINES][0]
                    point = get_intersection_point(other_line[0], other_line[1], -1, 1, precision=15)
                else:
                    point = get_intersection_point(lines_comb[0][0], lines_comb[0][1], lines_comb[1][0],
                                                   lines_comb[1][1], precision=15)
                if point:
                    if point in nodes:
                        nodes[point].update(lines_comb)
                    else:
                        nodes[point] = set(lines_comb)
    return nodes


def merge_nodes_geometrically_recursion(nodes, merge, tag):
    """
    Merge nodes recursively according to geometrical restrictions
    """
    triangle_lines = ('x_axis' or 'hypotenuse' or 'y_axis')
    for n1, n2 in combinations(sorted(nodes), 2):
        if len(nodes[n1] & nodes[n2]) >= 2:
            tag = True
            merged = nodes[n1].union(nodes[n2])
            if (triangle_lines in nodes[n2] and triangle_lines not in nodes[n1]) or len(nodes[n1]) < len(nodes[n2]):
                nodes[n2] = merged
                nodes[n1] = set()
                merge[n1] = n2
            else:
                nodes[n1] = merged
                nodes[n2] = set()
                merge[n2] = n1
            return merge_nodes_geometrically_recursion(nodes, merge, tag)
    return nodes, merge, tag


def merge_nodes_geometrically(nodes, tag):
    """
    Merge nodes according to geometrical restrictions
    """
    geometry_nodes, merge, tag = merge_nodes_geometrically_recursion(nodes, {}, tag)
    geometry_nodes = {k: l for k, l in geometry_nodes.items() if len(l) != 0}
    return geometry_nodes, merge, tag


def get_merged_areas_recursion(merge, areas):
    """
    Merge nodes recursively according to area superposition
    """
    for index, n1 in enumerate(merge):
        if len(merge) != 1:
            n2 = merge[index + 1]
            if min(n2[0]) <= max(n1[0]):
                if min(n1[1]) <= min(n2[1]) <= max(n1[1]) or min(n1[1]) <= max(n2[1]) <= max(n1[1]):
                    merge.insert(index, ((min(n1[0] + n2[0]), max(n1[0] + n2[0])),
                                         (min(n1[1] + n2[1]), max(n1[1] + n2[1]))))
                    merge.pop(index + 1)
                    merge.pop(index + 1)
                    return get_merged_areas_recursion(merge, areas)
        areas.append(n1)
        merge.pop(0)
        if merge:
            return get_merged_areas_recursion(merge, areas)
    return areas


def get_merged_areas(merge):
    """
    Get the area of the nodes merged geometrically
    """
    x_values = []
    y_values = []
    for old, new in merge.items():
        old_x, old_y = old
        new_x, new_y = new
        x_values.append((old_x, new_x))
        y_values.append((old_y, new_y))
    merged_nodes_sorted = sorted(zip(x_values, y_values), key=lambda x: min(x[0]))
    areas = get_merged_areas_recursion(merged_nodes_sorted, [])
    return areas


def get_nodes_in_area_recursion(area_nodes, geometry_nodes, nodes_analyse, area):
    """
    Get the node corresponding to the merged area
    """
    if len(nodes_analyse) == 1:
        return geometry_nodes
    for n1, n2 in combinations(nodes_analyse, 2):
        merged = geometry_nodes[n1].union(geometry_nodes[n2])
        if len(TRIANGLE_LINES.union(geometry_nodes[n1])) < len(TRIANGLE_LINES.union(geometry_nodes[n2])) \
                or len(geometry_nodes[n1]) < len(geometry_nodes[n2]):
            area_nodes[area].remove(n1)
            nodes_analyse.remove(n1)
            geometry_nodes[n2] = merged
            geometry_nodes[n1] = set()
        else:
            area_nodes[area].remove(n2)
            nodes_analyse.remove(n2)
            geometry_nodes[n1] = merged
            geometry_nodes[n2] = set()
        get_nodes_in_area_recursion(area_nodes, geometry_nodes, nodes_analyse, area)


def get_nodes_in_area(areas, merge, geometry_nodes):
    """
    Get the node corresponding to the merged area
    """
    new_nodes = merge.values() - merge.keys()
    area_nodes = defaultdict(set)
    for area in areas:
        for node in new_nodes:
            x, y = node
            if min(area[0]) <= x <= max(area[0]) and min(area[1]) <= y <= max(area[1]):
                area_nodes[area].add(node)
    for area, nodes_to_analyse in area_nodes.items():
        get_nodes_in_area_recursion(area_nodes, geometry_nodes, nodes_to_analyse, area)
    return area_nodes


def merge_nodes_in_area(area_nodes, geometry_nodes):
    """
    Check in the areas if there are other nodes that have been previously merged
    """
    merged = {}
    for area, real_node in area_nodes.items():
        real_node = real_node.pop()
        if real_node not in merged:
            merged[real_node] = geometry_nodes[real_node]
        for node, lines in geometry_nodes.items():
            x, y = node
            if min(area[0]) <= x <= max(area[0]) and min(area[1]) <= y <= max(area[1]):
                merged[real_node].update(geometry_nodes[real_node].union(lines))
            else:
                if node not in merged.keys():
                    merged[node] = lines
    return merged


def merge_nodes(nodes):
    """
    Merge nodes by geometrical rules
    """
    tag = False
    geometry_nodes, merge, tag = merge_nodes_geometrically(nodes, tag)
    if tag:
        areas = get_merged_areas(merge)
        area_nodes = get_nodes_in_area(areas, merge, geometry_nodes)
        merged_nodes = merge_nodes_in_area(area_nodes, geometry_nodes)
        return merge_nodes(merged_nodes)
    return geometry_nodes


def get_lines(nodes):
    """
    Get lines according to merged nodes
    """
    lines = defaultdict(list)
    for node, node_lines in nodes.items():
        for line in node_lines:
            lines[line].append(node)
    return lines


def sort_line_nodes(lines):
    """
    Sort nodes according to their distance from the y or x axis
    """
    sorted_lines = {}
    for line, nodes in lines.items():
        temporal_x = [n[0] for n in nodes]
        temporal_y = [n[1] for n in nodes]
        if line == 'y_axis':
            origin = min(nodes, key=lambda n: n[1])
        else:
            if temporal_x.count(0.0) == 1 or temporal_y.count(0.0) == 0:
                origin = min(nodes, key=lambda n: n[0])
            elif temporal_y.count(0.0) == 1 or temporal_x.count(0.0) == 0:
                origin = min(nodes, key=lambda n: n[1])
            else:
                raise Exception('ERROR: nodes of line ' + str(line) + " couldn't be sorted with nodes " + str(nodes))
        distances_temp = {node: math.hypot(node[0] - origin[0], node[1] - origin[1]) for node in nodes}
        sorted_lines[line] = [n for n, distance in sorted(distances_temp.items(), key=lambda x: x[1])]
    return sorted_lines


def get_predictors_intersection(rho, predictors, precision):
    """
    Get lines and nodes of predictors's planes intersection
    """
    nodes = get_nodes(rho, predictors, precision)
    nodes = unmerge_nodes(nodes)
    nodes = merge_nodes(nodes)
    lines = get_lines(nodes)
    sorted_lines = sort_line_nodes(lines)
    return nodes, sorted_lines
