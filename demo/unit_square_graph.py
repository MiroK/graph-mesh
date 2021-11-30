from graph_mesh.generators import random_unit_square
from graph_mesh import *
from dolfin import File
from xii.meshing.embedded_mesh import TangentCurve
    
G, edge_f = random_unit_square(32, prob=0.7)
mesh, radii_f = mesh_graph(G)

File('edge_f.pvd') << edge_f
File('radii_f_coarse.pvd') << radii_f    
# Refine the mesh
for _ in range(4):
    rmesh = graph_adapt(mesh)
    rradii_f = graph_adapt(radii_f, rmesh)
        
    mesh, radii_f = rmesh, rradii_f

File('radii_f_fine.pvd') << radii_f

cell_f, bcolors, lcolors, terminals = color_branches(mesh)

File('graph_colored.pvd') << cell_f

tangent = TangentCurve(mesh)
File('tangent.pvd') << tangent
