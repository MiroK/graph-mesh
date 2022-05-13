def trim_soma(graph):
    '''Remove the soma'''
    # Soma of neuron is often represented as cylinder which in the centerline
    # representation has more segments. We want to reduce it to a point
    remove = [edge for edge in graph.edges if graph.get_edge_data(*edge)['ntype'] == 1]
    remove and [graph.remove_edge(*e) for e in remove]

    return graph

# ----------------------------------------------------------------------------

if __name__ == '__main__':
    from graph_mesh.embedding import box_embed
    from graph_mesh.swc import swc2graph
    import dolfin as df
    import numpy as np
    import sys
    
    swc_path = 'RatS1-10-91.CNG.swc'
    graph = swc2graph(swc_path)
    graph = trim_soma(graph)

    for i in range(5):
        scale = 0.4/2**i
        args = ['', '-clscale', str(scale)]
        mesh, entity_fs = box_embed(graph, scaling=np.array([1.1, 1.1, 1.1]), args=args)

        with df.HDF5File(mesh.mpi_comm(), f'neuron_{i}.h5', 'w') as out:
            out.write(mesh, '/mesh')
            out.write(entity_fs[1], '/neuron')
