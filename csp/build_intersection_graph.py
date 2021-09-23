"""
Get the node's interactions, expected edges and expected nodes in all of polygons
"""

from collections import defaultdict
from find_predictor_intersections import nodes2edge


def get_search_edges(lines):
    """
    Get the segments/edges between nodes of a line
    """
    search_edges = []
    edge_type = {}
    for line, nodes in lines.items():
        exterior_line = False
        if line in ['x_axis', 'y_axis', 'hypotenuse']:
            exterior_line = True
        for index, node in enumerate(nodes[:-1]):
            edge = nodes2edge(node, nodes[index + 1])
            if exterior_line:
                search_edges.append(edge)
                edge_type[edge] = 'exterior'
            else:
                search_edges.extend([edge] * 2)
                edge_type[edge] = 'interior'
    return sorted(search_edges), edge_type


def get_interactions(search_edges):
    """
    Get the interaction between nodes
    """
    interactions = defaultdict(list)
    for edge in sorted(set(search_edges)):
        n1, n2 = edge
        interactions[n1].append(n2)
        interactions[n2].append(n1)
    return interactions


def get_search_nodes(interactions, edge_type):
    """
    Get the number of polygons that we expect for each node
    that is, in how many polygons a node will appear
    """
    search_nodes = []
    for n1, children in interactions.items():
        expected_polygons = len(children)
        for n2 in children:
            edge = nodes2edge(n1, n2)
            if edge_type[edge] == 'exterior':
                expected_polygons -= 1
                break
        search_nodes.extend([n1] * expected_polygons)
    return sorted(search_nodes)


def get_predictors_graph(lines):
    """
    Get the predictors graph
    """
    search_edges, edge_type = get_search_edges(lines)
    interactions = get_interactions(search_edges)
    search_nodes = get_search_nodes(interactions, edge_type)
    return interactions, search_edges, search_nodes
