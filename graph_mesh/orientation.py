from xii import EmbeddedMesh
from xii.meshing.embedded_mesh import TangentCurve
import numpy as np
from dolfin import *
from xii import *


def compute_io_orientation(color_f):
    '''
    Let color_f be a cell function of a 1d in (dim > 1) mesh marking the 
    branches. Using a mesh global tangent field tau we compute facet functions
    of the branch such that for facet marked with 1 there is dot(tau, n) > 0 
    and 2 otherwise. 

    Return {color -> (oriented facet function (vertex functio), P1)}
    '''
    # NOTE: computed tangent orientation does not necessarily coincide with
    # fenics only so fenics tangential derivative should probably be replaced by
    # dot(grad(f), tau)
    mesh = color_f.mesh()

    tau = TangentCurve(mesh)
    tangents = tau.vector().get_local()

    # What we want to build
    oriented = {}

    branch_colors = np.unique(color_f.array())
    for color in branch_colors:
        branch = EmbeddedMesh(color_f, color)  # This should be MeshView in FEniCS?
        dx_ = Measure('dx', domain=branch)
        
        # Localize tangent
        V = VectorFunctionSpace(branch, 'DG', 0)
        v = TestFunction(V)

        vol = CellVolume(branch)
        
        tau_branch = Function(V)
        tau_branch.vector()[:] = ii_assemble((1/vol)*inner(Restriction(tau, branch), v)*dx_)

        n = FacetNormal(branch)
        # Indicator space for each branch
        V = FunctionSpace(branch, 'CG', 1)
        v = TestFunction(V)

        area = FacetArea(branch)
        f = Function(V)
        # Test conormality with our tangent
        io_condition = conditional(gt(dot(tau_branch, n), Constant(0)), Constant(1), Constant(2))
        o = assemble((1/area)*inner(v, io_condition)*ds)
        f.vector()[:] = o
        
        o = np.asarray(o, int)
        # Now it remains to make P1 into a facet function
        facet_f = MeshFunction('size_t', branch, 0, 0)
        facet_f.array()[:] = o[vertex_to_dof_map(V)]
        
        oriented[color] = (facet_f, f)  # P1 is here just for plotting

    return oriented
