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
    from graph_mesh.coloring import greedy_color
    from graph_mesh.swc import swc2graph
    import dolfin as df
    import numpy as np
    

    subject = 'BG0002'
    swc_path = get_swc(subject)
    G = swc2graph(swc_path)

    mesh, radii_f, ntype_f = mesh_graph(G)

    File(f'{subject}_radii_f_coarse.pvd') << radii_f
    
    cell_f, bcolors, lcolors, terminals = color_branches(mesh)
    File(f'{subject}_graph_colored.pvd') << cell_f

    sparse_color = greedy_color(cell_f, terminals)
    File(f'{subject}_sparse_color.pvd') << sparse_color

    # Save mesh as coordinates and cell 2 vertex matrix
    coordinates = mesh.coordinates()
    v2c = mesh.cells()    # Each row is a cell encoded in terms of its vertices

    from scipy.io import savemat, loadmat

    savemat('f{subject}.mat', {'coordinates': coordinates, 'cells': v2c})

    data = loadmat('f{subject}.mat')

    assert np.linalg.norm(coordinates-data['coordinates']) < 1E-13
    assert np.linalg.norm(v2c-data['cells']) < 1E-13    

    # Just for illustration
    from graph_mesh.orientation import compute_io_orientation
    from xii import *
    import time
    
    tau = TangentCurve(mesh)

    submesh_data = compute_io_orientation(cell_f, tau)
    
    # There will be one global pressure space
    Q = df.FunctionSpace(mesh, 'CG', 1)
    p, q = df.TrialFunction(Q), df.TestFunction(Q)

    Grad = lambda f, t: df.dot(df.grad(f), t)

    A01, A10 = 0, 0

    print('PSEUDOSTOKES ASSEMBLY')
    then = time.time()

    nbranches = len(bcolors)
    # This is withough the coupling
    a = np.zeros((nbranches+1, nbranches+1)).tolist()

    Qindex = nbranches
    
    for Vindex, color in enumerate(bcolors):
        facet_f, branch, tau_branch = submesh_data[color]

        Vi = df.FunctionSpace(branch, 'CG', 2)
        ui, vi = df.TrialFunction(Vi), df.TestFunction(Vi)

        Rp_i, Rq_i = Restriction(p, branch), Restriction(q, branch)
        dxi = df.Measure('dx', domain=branch)

        a[Vindex][Vindex] = df.inner(Grad(ui, tau_branch), Grad(vi, tau_branch))*dxi
        a[Vindex][Qindex] = df.inner(Grad(vi, tau_branch), Rp_i)*dxi
        a[Qindex][Vindex] = df.inner(Grad(ui, tau_branch), Rq_i)*dxi

    A = ii_assemble(a)

    print('Done', time.time()-then)

    for Vindex, color in enumerate(bcolors):
        assert A[Vindex][Vindex].norm('linf') > 0
        assert ii_convert(A[Vindex][Qindex]).norm('linf') > 0
        assert ii_convert(A[Qindex][Vindex]).norm('linf') > 0        
