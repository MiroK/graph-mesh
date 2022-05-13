from gmshnics import msh_gmsh_model, mesh_from_gmsh
import numpy as np
import gmsh


def get_bbox(graph, scaling, align):
    '''Bounding box'''
    assert np.all(scaling >= 1)
    
    nodes = graph.nodes
    mesh_coordinates = np.row_stack([nodes[n]['pos'] for n in nodes])

    if align == False:
        ll, uu = np.min(mesh_coordinates, axis=0), np.max(mesh_coordinates, axis=0)

        dx = uu - ll
        shift = 0.5*(scaling - 1)
        origin = ll - shift*dx
        dx = scaling*dx

        return (origin,
                dx[0]*np.array([1, 0, 0]),
                dx[1]*np.array([0, 1, 0]),
                dx[2]*np.array([0, 0, 1]))

    raise ValueError


def box_embed(graph, scaling, align=False, args=[]):

    origin, dx, dy, dz = get_bbox(graph, scaling, align=align)

    gmsh.initialize(args)

    model = gmsh.model
    fac = model.geo
    # Start from the outside
    entities = addBox(fac, origin, dx, dy, dz)

    nodes = graph.nodes
    vertices = [fac.addPoint(*nodes[n]['pos']) for n in nodes]

    lines = [fac.addLine(vertices[p-1], vertices[q-1]) for p, q in graph.edges]
    fac.synchronize()
    # Mark them for lookup
    model.addPhysicalGroup(1, lines, 1)
    # Mark volume to get it meshes too
    vol = entities[3][-1]
    model.addPhysicalGroup(3, [vol], 1)
    # Embedd the lines
    model.mesh.embed(1, lines, 3, vol)    
    fac.synchronize()

    nodes, topologies = msh_gmsh_model(model, 3)
    mesh, entity_functions = mesh_from_gmsh(nodes, topologies)    

    gmsh.finalize()

    return mesh, entity_functions


def addBox(fac, origin, dx, dy, dz):
    '''This is
     7-----6
     /|   /|  
    4----5 |
    | 3--|-2
    |/   |/     /
    0----1      --->x
    '''
    points = [origin,
              origin + dx,
              origin + dx + dy,
              origin + dy,
              origin + dz,
              origin + dz + dx,
              origin + dz + dx + dy,
              origin + dz + dy]
    points = [fac.addPoint(*p) for p in points]

    lines = [(0, 1), (1, 2), (2, 3), (3, 0),
             (0, 4), (1, 5), (2, 6), (3, 7),
             (4, 5), (5, 6), (6, 7), (7, 4)]

    lines = [fac.addLine(points[p], points[q]) for p, q in lines]

    surfaces = [(7, 11, -4, -3), (1, 6, -9, -5),
                (0, 5, -8, -4), (6, 10, -7, -2),
                (0, 1, 2, 3), (8, 9, 10, 11)]

    plane_surfaces = []
    for surf in surfaces:
        cl = fac.addCurveLoop([lines[s] if s >= 0 else -lines[-s] for s in surf])
        plane_surfaces.append(fac.addPlaneSurface([cl]))
    sl = fac.addSurfaceLoop(plane_surfaces)
    vol = fac.addVolume([sl])
    
    fac.synchronize()

    return {3: [vol], 2: plane_surfaces}


# --------------------------------------------------------------------


if __name__ == '__main__':
    import numpy as np
    import sys, gmsh

    gmsh.initialize(sys.argv)

    model = gmsh.model
    fac = model.geo

    origin = np.zeros(3)
    dx, dy, dz = np.eye(3)
    
    entities = addBox(fac, origin, dx, dy, dz)

    center = fac.addPoint(0.5, 0.5, 0.5)
    top = fac.addPoint(0.5, 0.5, 1)
    
    l = fac.addLine(center, top)

    fac.synchronize()
    
    model.mesh.embed(0, [top], 2, entities[2][-1])
    model.mesh.embed(1, [l], 3, entities[3][-1])        

    gmsh.fltk.initialize()
    gmsh.fltk.run()

    gmsh.finalize()
