# Here we work with SWC files encoding reconstructions of networks from
# BRAVA DATABASE http://cng.gmu.edu/brava/home.php?s=1&name_browser=false

# The Brain Vasculature (BraVa) database contains digital reconstructions of
# the human brain arterial arborizations from 61 healthy adult subjects along
# with extracted morphological measurements. The arterial arborizations include
# the six major trees stemming from the circle of Willis, namely: the left and
# right Anterior Cerebral Arteries (ACAs), Middle Cerebral Arteries (MCAs), and
# Posterior Cerebral Arteries (PCAs).
import os
import wget

def get_swc(subject):
    '''Download it'''
    assert not subject.endswith('.swc')

    swc_path = f'{subject}.CNG.swc'
    if os.path.exists(swc_path):
        return swc_path

    url = f'http://cng.gmu.edu/brava/files/file_swc/original/{swc_path}'
    swc_path = wget.download(url)

    return swc_path

# --------------------------------------------------------------------

if __name__ == '__main__':
    from graph_mesh import *
    from dolfin import File
    from xii.meshing.embedded_mesh import TangentCurve
    from graph_mesh.orientation import compute_io_orientation
    from graph_mesh.coloring import greedy_color
    from graph_mesh.swc import swc2graph

    subject = 'BG001'
    swc_path = get_swc(subject)
    G = swc2graph(swc_path)

    mesh, radii_f = mesh_graph(G)

    File(f'{subject}_radii_f_coarse.pvd') << radii_f
    
    cell_f, bcolors, lcolors, terminals = color_branches(mesh)
    File(f'{subject}_graph_colored.pvd') << cell_f

    tangent = TangentCurve(mesh)
    File(f'{subject}_tangent.pvd') << tangent

    sparse_color = greedy_color(cell_f, terminals)
    File(f'{subject}_sparse_color.pvd') << sparse_color