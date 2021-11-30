import networkx as nx
import numpy as np
import networkgen.geometrytools  as geometrytools
import networkgen.generationtools as generationtools
import random


# https://gitlab.com/ValletAlexandra/NetworkGen/-/blob/master/scripts/Tree.py
def make_graph(N):
    '''Creation of an tree-like network'''
    # Parameters
    # Origin location
    p0=[0,0,0]

    # Initial direction
    direction=[0,1,0]

    # First vessel diameter
    D0=1

    # lambda
    lmbda=8

    # gamma
    gam=0.8

    # By convention we chose gam <=1 so D1 will always be smaller or equal to D2
    if gam > 1:
        raise Exception('Please choose a value for gamma lower or equal to 1')


    #Surface normal function
    #The surface normal here is fixed because we want to stay in the x,y plane. 
    #But this could be the normal of any surface.
    def normal (x,y,z) :
        return [0,0,1]

    #### Creation of the graph

    # Create a networkx graph 
    G=nx.DiGraph()

    # Create the first vessel
    #########################
    L=D0*lmbda
    print('Root :')
    print('D :',D0)
    print('L:',L)


    G.add_edge(0,1)

    nx.set_node_attributes(G, p0, "pos")
    nx.set_edge_attributes(G, D0/2, "radius")


    G.nodes[1]['pos']=geometrytools.translation(p0,direction,L)
    inode=1

    #### Iteration to create the other vessels following a given law

    #list of the vessels from the previous generation
    previous_edges=[(0,1)]

    for igen in range(1,N):
        current_edges=[]
        for e in previous_edges : 
            # Parent vessel properties
            previousvessel=[G.nodes[e[0]]['pos'],G.nodes[e[1]]['pos']]
            D0=G.edges[e]['radius']*2

            # Daughter diameters
            D2=D0*(gam**3+1)**(-1/3)
            D1=gam*D2
            # Daughter lenghts
            L2=lmbda*D2
            L1=lmbda*D1
            # Bifurcation angles
            # angle for the smallest vessel
            cos1= (D0**4 +D1**4 -(D0**3 - D1**3)**(4/3))/(2*D0**2*D1**2)
            angle1= np.degrees(np.arccos(cos1))
            # angle for the biggest vessel
            cos2=(D0**4 +D2**4 -(D0**3 - D2**3)**(4/3))/(2*D0**2*D2**2)
            angle2=np.degrees(np.arccos(cos2))

            #randomly chose which vessel go to the right/left
            sign1=random.choice((-1, 1))
            sign2=-1*sign1


            ### Add first daughter vessel
            print('Daughter 1 :')
            print('D :',D1)
            print('L:',L1)
            print('angle :',angle1)

            inode+=1
            new_edge=(e[1],inode)
            G.add_edge(*new_edge)

            # Set the location according to length and angle
            G.nodes[inode]['pos']=generationtools.compute_vessel_endpoint (previousvessel, normal(*previousvessel[1]),sign1*angle1,L1)

            # Set radius
            G.edges[new_edge]['radius']=D1/2

            # Add to the pool of vessels for this generation
            current_edges.append(new_edge)

            ### Add second daughter vessel
            print('Daughter 2 :')
            print('D :',D2)
            print('L :',L2)
            print('angle :',angle2)

            inode+=1
            new_edge=(e[1],inode)
            G.add_edge(*new_edge)

            # Set the location according to length and angle
            G.nodes[inode]['pos']=generationtools.compute_vessel_endpoint (previousvessel, normal(*previousvessel[1]),sign2*angle2,L2)

            # Set radius
            G.edges[new_edge]['radius']=D2/2

            # Add to the pool of vessels for this generation
            current_edges.append(new_edge)
        previous_edges=current_edges

    return G

# ----------------------------------------------------------------------------

if __name__ == '__main__':
    from graph_mesh import *
    from dolfin import File
    from xii.meshing.embedded_mesh import TangentCurve
    from graph_mesh.orientation import compute_io_orientation
    
    # NOTE: one issue with this mesh generator is that it allows children
    # branches to cross parants. However, this happens if there are many generations
    # so here we keep it low to be on the safe side
    G = make_graph(N=5)
    mesh, radii_f = mesh_graph(G)

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

    # Branches with IO
    oriented = compute_io_orientation(cell_f)

    for color, (ff, foo) in oriented.items():
        File(f'ff_{color}.pvd') << foo
