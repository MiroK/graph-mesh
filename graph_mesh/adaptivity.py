import dolfin as df


def graph_adapt(first, second=None):
    '''Extend adapt to hand DG0 functions'''
    if second is None:
        return df.adapt(first)

    try:
        return df.adapt(first, second)
    except TypeError:
        assert first.function_space().ufl_element().family() == 'Discontinuous Lagrange'
        # Second assumed to be refined mesh of first
        rmesh = second
        
        V = first.function_space()
        dm = V.dofmap()
        values = first.vector().get_local()

        tdim = rmesh.topology().dim()
        mapping = rmesh.data().array('parent_cell', tdim)

        rV = df.FunctionSpace(rmesh, V.ufl_element())
        rdm = rV.dofmap()
        # This we want to fill
        f = df.Function(rV)
        rvalues = f.vector().get_local()
        # Set my values taking those from the parent
        for rc, c in enumerate(mapping):
            rvalues[rdm.cell_dofs(rc)] = values[dm.cell_dofs(c)]

        f.vector().set_local(rvalues)

        return f
