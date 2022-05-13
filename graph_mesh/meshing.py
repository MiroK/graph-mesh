from dolfin import MeshEditor, Mesh, Function, FunctionSpace
import networkx as nx
import numpy as np


def mesh_graph(graph):
    '''1d mesh of graph. Assume nodes have data `pos` and edges have `radius`'''
    nodes = graph.nodes

    mesh_coordinates = np.row_stack([nodes[n]['pos'] for n in nodes])
    nvtx, gdim = mesh_coordinates.shape
    tdim = 1

    mapping = dict((n, i) for i, n in enumerate(nodes))
    
    mesh = Mesh()
    editor = MeshEditor()
    editor.open(mesh, 'interval', tdim, gdim)            

    editor.init_vertices(nvtx)
    editor.init_cells(graph.number_of_edges())

    for i, x in enumerate(mesh_coordinates):
        editor.add_vertex(i, x)

    radius_cell_data = np.zeros(graph.number_of_edges())
    ntype_cell_data = np.zeros(graph.number_of_edges())

    edges = graph.edges()
    for ci, (ni, nj) in enumerate(edges):
        editor.add_cell(ci, (mapping[ni], mapping[nj]))
        radius_cell_data[ci] = edges[(ni, nj)]['radius']
        ntype_cell_data[ci] = edges[(ni, nj)]['ntype']        
    
    editor.close()

    radius_data = Function(FunctionSpace(mesh, 'DG', 0))
    radius_data.vector().set_local(radius_cell_data)

    ntype_data = Function(FunctionSpace(mesh, 'DG', 0))
    ntype_data.vector().set_local(ntype_cell_data)

    return mesh, radius_data, ntype_data
