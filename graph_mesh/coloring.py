from collections import defaultdict
import networkx as nx
import numpy as np
import dolfin as df


def greedy_color(color_f, color_connectivity):
    '''
    Let cell_f mark uniquely the mesh branches. As this might be many 
    colors the coloring computed here strives to reduce that by having 
    no two branches that are connected be colored with the same color
    '''
    mesh = color_f.mesh()
    assert mesh.geometry().dim() > 1 and mesh.topology().dim() == 1
    assert color_f.dim() == mesh.topology().dim()
    
    dG = nx.Graph()
    dG.add_edges_from(color_connectivity[c] for c in sorted(color_connectivity))
    # We want to build a graph where branch is a node and edges represent
    # branch connectivity.
    G = nx.line_graph(dG)
    # Then the job is to color the graph ...
    colormap = nx.algorithms.coloring.greedy_color(G)
    
    # Unique coloring from branch encoding in terms of vertices
    icc = {tuple(sorted(v)): c for c, v in color_connectivity.items()}

    color_f = color_f.array()
    sparse_color = 1*color_f  # A copy

    greedy_colors = set()
    # ... and translate the coloring
    for key, color in colormap.items():
        branch, = np.where(color_f == icc[key])
        sparse_color[branch] = color+1
        # Just to keep track of how many we ended up using
        greedy_colors.add(color+1)

        del icc[key]
    assert not icc

    ans = df.MeshFunction('size_t', mesh, 1, 0)
    ans.array()[:] = sparse_color

    print(f'Greedy brach coloring color reduction {len(color_connectivity)} -> {len(greedy_colors)}')
    
    return ans
