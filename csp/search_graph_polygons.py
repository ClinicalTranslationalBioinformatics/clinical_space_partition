"""
Get the polygons in which the the cost space is divided by the predictors
"""

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from find_predictor_intersections import nodes2edge
from itertools import combinations


def get_polygon(first_node, nodes, interactions, search_edges, paths, covered_lines, found_polygons):
    """
    Get a polygon
    """
    path = paths.pop(0)
    node = path[-1]
    path_lines = covered_lines.pop(0)
    node_edge_not_in_polygon = True
    if path_lines and node == first_node:
        for every_node in nodes:
            if every_node not in path:
                # C6: Discard because there is a point inside the polygon
                if Polygon(path).contains(Point(every_node)):
                    every_node_interactions_polygon_nodes = [every_node in interactions[n] for n in path[:-1]]
                    if every_node_interactions_polygon_nodes.count(True) >= 4:
                        node_edge_not_in_polygon = False
                        break
                else:
                    # C7: Discard because there is an edge inside the polygon
                    polygon_edges = [nodes2edge(n, path[i + 1]) for i, n in enumerate(path[:-1])]
                    for n1, n2 in combinations(path[:-1], 2):
                        edge = nodes2edge(n1, n2)
                        if edge in search_edges and edge not in polygon_edges:
                            node_edge_not_in_polygon = False
                            break
    if path_lines and node == first_node and node_edge_not_in_polygon:
        return path
    else:
        for child in interactions[node]:
            current_line = nodes[node] & nodes[child]
            edge = nodes2edge(node, child)
            # C1: Discard because the edge is in a line already seen
            if len(current_line & path_lines) > 0:
                continue
            # C2: Discard because the node is already in the path and is not the first node
            elif child in path and child != first_node:
                continue
            # C3: Discard because the edge can't participate in more polygons
            elif edge not in search_edges:
                continue
            # C5: Discard because the polygon has already been found
            elif set(path + [child]) in found_polygons:
                continue
            covered_lines.append(path_lines.union(current_line))
            paths.append(path + [child])
        return get_polygon(first_node, nodes, interactions, search_edges, paths, covered_lines, found_polygons)


def get_polygons(search_edges, search_nodes, nodes, interactions):
    """
    Get the polygons
    """
    polygons = []
    found_polygons = []
    while search_nodes:
        first_node = search_nodes[0]
        polygon = get_polygon(first_node, nodes, interactions, search_edges, [[first_node]], [set()], found_polygons)
        polygons.append(polygon)
        found_polygons.append(set(polygon))
        # C4: Decrease polygon's nodes counters
        for node in polygon[:-1]:
            search_nodes.remove(node)
        # C3: Decrease polygon's edges counters
        for index, node in enumerate(polygon[:-1]):
            edge = nodes2edge(node, polygon[index + 1])
            search_edges.remove(edge)
    return polygons
