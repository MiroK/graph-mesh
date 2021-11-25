import networkx as nx
import dolfin as df
import numpy as np


def random_edge_graph(mesh, prob=0.5):
    '''Graph based on the edge function'''
    assert mesh.topology().dim() > 1
    
    edge_f = df.MeshFunction('size_t', mesh, 1, 0)
    num_edges = mesh.num_entities(1)
    
    edge_f.array()[np.random.rand(num_edges) > prob] = 1

    mesh.init(1, 0)
    
    G = nx.Graph()
    radii = np.random.rand(sum(edge_f.array() == 1))
    G.add_edges_from(tuple(e.entities(0)) + ({'radius': r}, )
                      for e, r in zip(df.SubsetIterator(edge_f, 1), radii))
    # Attach positions
    x = mesh.coordinates()

    nodes = G.nodes
    for i in nodes:
        nodes[i]['pos'] = x[i]

    return G, edge_f


def random_unit_square(N, prob=0.5):
    '''Based on (0, 1)^2'''
    mesh = df.UnitSquareMesh(N, N)

    return random_edge_graph(mesh, prob=prob)


