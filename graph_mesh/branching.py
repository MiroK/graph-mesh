from copy import deepcopy
import dolfin as df


def color_branches(mesh):
    '''Start/end is a terminal

    We return a cell function where branches and loops are colored,
    a list of brach colors, list of loop colors and map from terminal 
    vertices to connected cells.
    '''
    assert mesh.topology().dim() == 1

    _, v2c = mesh.init(0, 1), mesh.topology()(0, 1)
    c2v = mesh.topology()(1, 0)

    terminals = {v: set(v2c(v)) for v in range(mesh.num_vertices()) if len(v2c(v)) != 2}
    
    cell_f = df.MeshFunction('size_t', mesh, 1, 0)
    if not terminals:
        cell_f.set_all(1)
        return cell_f, [], [1]

    def next_vertex(c, v, c2v=c2v):
        v0, v1 = c2v(c)
        return v1 if v == v0 else v0

    def next_cell(v, c, v2c=v2c):
        c0, c1 = v2c(v)
        return c1 if c == c0 else c0
    
    values = cell_f.array()
    branch_colors, loop_colors, color = [], [], 0
    # It will be useful to talk about end vertices of colored branches
    color_connectivity = {}

    exhausted = False
    while not exhausted:
        vertex = max(terminals, key=lambda v: terminals[v])
        vertex_cells = terminals[vertex]

        exhausted = len(vertex_cells) == 0
        # The idea is to walk from vertex following the cell
        while vertex_cells:
            link_cell = vertex_cells.pop()
            v0 = vertex
            
            branch = [link_cell]
            # v0 --
            while next_vertex(link_cell, v0) not in terminals:
                # -- v0 ==
                v0 = next_vertex(link_cell, v0)
                # Because we have not terminal, ==
                link_cell = next_cell(v0, link_cell)
                branch.append(link_cell)
            # Once we reached the terminal
            v0 = next_vertex(link_cell, v0)

            color += 1
            if v0 == vertex:
                loop_colors.append(color)
            else:
                branch_colors.append(color)
            values[branch] = color

            color_connectivity[color] = (v0, vertex)
            
            # Preclude leaving from vertex in a look
            link_cell in vertex_cells and vertex_cells.remove(link_cell)
            # If we arrived to some other terminal, we don't want to leave from it by the
            # same way we arrived
            v0 in terminals and link_cell in terminals[v0] and terminals[v0].remove(link_cell)

    return cell_f, branch_colors, loop_colors, color_connectivity
