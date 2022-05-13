import networkx as nx
import numpy as np
import itertools


def swc2graph(swc_file):
    '''Parse SWC file to valid graph for `graph-mesh`'''
    parse = lambda l: [fi(li) for fi, li in zip((int, int, float, float, float, float, int), l)]
    # 1 1 0.0 0.0 0.0 7.3875 -1

    with open(swc_file, 'r') as swc:
        comments = []
        # Get comments out of the way
        iter_points = map(lambda l: parse(l.strip().split()),
                          # Talk about lambdas without side-effects :)
                          itertools.dropwhile(lambda l: l.startswith('#') and comments.append(l) is None, swc))
        # This is the first point - can't make pairs yet
        index, ntype, x, y, z, radius, parent = next(iter_points)
        # SWC starts numbering at 1
        assert index == 1

        G = nx.Graph()
        G.add_node(index, pos=np.array([x, y, z]))
                   
        for count, (index, ntype, x, y, z, radius, parent) in enumerate(iter_points, index+1):

            assert count == index

            G.add_node(index, pos=np.array([x, y, z]))
            G.add_edge(index, parent, radius=radius, ntype=ntype)

    return G
