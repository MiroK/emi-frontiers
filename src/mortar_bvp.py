# Here we solve the coupled Laplacian problem (*) arising in the time
# discretization of the mortar formulation of the EMI model by backward
# Euler scheme:
#
#     -div(sigma_i nabla u_i) = f_i on (0.25, 0.75)^2,
#     -div(sigma_e nabla u_e) = f_e on (0, 1)^2\(0.25, 0.75)^2
#
# With Dirichlet boundary conditions
#
#     u_e = g_e on for each of bdries {1, 2, 3, 4}
#
# And coupling conditions on bdries 
# 
#     dot(ni, sigma_i nabla u_i) + dot(ne, sigma_e nabla u_e) = 0,
#     u_i - u_e -beta*p = f_m
#
# here p = -dot(ni, sigma_i nabla u_i)
# 
#         3
#  -------------------
#  |     13          |
#  |    |--------|   |
# 4| 14 |  i(10) |12 |  2
#  |    |--------|   |
#  |e(1)  11         |
#  -------------------
#        1
import numpy as np
from block import block_bc, block_vec
from dolfin import *
from xii import *


def has_sane_inputs(mesh, sigma, f, beta, bcs):
    '''Check sanity of inputs'''
    # We have [0, 1]^2
    assert np.linalg.norm(mesh.coordinates().min(axis=0) - np.zeros(2)) < 1E-13
    assert np.linalg.norm(mesh.coordinates().max(axis=0) - np.ones(2)) < 1E-13

    is_positive = lambda v: (isinstance(v, Constant) and v(0) > 0) or v > 0
    # Positive conductivities for both domains
    assert set(sigma.keys()) == set(('i', 'e'))
    assert all(map(is_positive, sigma.values()))

    # We have a forcing for each subdomain
    assert set(f.keys()) == set(('i', 'e', 'm'))

    # We have a boundary condition for each edge
    assert set(bcs.keys()) == set((1, 2, 3, 4))

    return True


def setup_domains(mesh):
    '''Interior, exterior, interfaces meshes and exterior facet function'''
    inside = CompiledSubDomain('(0.25-tol < x[0]) && (x[0] < 0.75+tol) && (0.25-tol < x[1]) && (x[1] < 0.75+tol)',
                               tol=1E-10)

    left_e = CompiledSubDomain('near(x[0], 0)')
    right_e = CompiledSubDomain('near(x[0], 1)')
    upper_e = CompiledSubDomain('near(x[1], 1)')
    lower_e = CompiledSubDomain('near(x[1], 0)')
    
    left_i = CompiledSubDomain('near(x[0], 0.25) && (0.25-tol < x[1]) && (x[1]< 0.75+tol)', tol=1E-10)
    right_i = CompiledSubDomain('near(x[0], 0.75) && (0.25-tol < x[1]) && (x[1]< 0.75+tol)', tol=1E-10)
    upper_i = CompiledSubDomain('near(x[1], 0.75) && (0.25-tol < x[0]) && (x[0]< 0.75+tol)', tol=1E-10)
    lower_i = CompiledSubDomain('near(x[1], 0.25) && (0.25-tol < x[0]) && (x[0]< 0.75+tol)', tol=1E-10)

    # Cell function for defining exterior and interior domains
    cell_f = MeshFunction('size_t', mesh, 2, 1)
    inside.mark(cell_f, 10)
    # Break it to meshes
    interior_mesh = SubMesh(mesh, cell_f, 10)
    exterior_mesh = SubMesh(mesh, cell_f, 1)

    # Facet function for marking the interface (viewed from inner)
    facet_f = MeshFunction('size_t', interior_mesh, 1, 0)

    for tag, subd in enumerate([lower_i, right_i, upper_i, left_i], 11):
        subd.mark(facet_f, tag)

    interface_mesh = EmbeddedMesh(facet_f, (11, 12, 13, 14))
    
    # Mark exterior for boundary conditions
    facet_f = MeshFunction('size_t', exterior_mesh, 1, 0)
    
    for tag, subd in enumerate([lower_e, right_e, upper_e, left_e], 1):
        subd.mark(facet_f, tag)

    return (interior_mesh, exterior_mesh, interface_mesh, facet_f)


class EmiBvpProblem(object):
    '''Get the linear system due to (*) and P1 elements.'''
    def __init__(self, mesh, sigma, rhs, beta, bcs, check_eigvals=False):    
        '''
        Compute system matrix
        
        Input
          mesh: [0, 1]^1
          sigma: dictionary 'i', 'e' to positive constants that are the conductivities
          rhs: dictionary 'i', 'e', 'm' to forces for interior, exterior, membrane
          beta: positive constant for the diagonal term in the multiplier equation
          bcs: dict 1, 2, 3, 4 to bdry values for the exterior boundaries
        '''
        assert has_sane_inputs(mesh, sigma, rhs, beta, bcs)

        interior_mesh, exterior_mesh, iface_mesh, exterior_facet_f = setup_domains(mesh)
    
        # Setup block bilinear and linear forms
        V1 = FunctionSpace(interior_mesh, 'CG', 1)
        V2 = FunctionSpace(exterior_mesh, 'CG', 1)
        Q = FunctionSpace(iface_mesh, 'CG', 1)
        W = [V1, V2, Q]

        u1, u2, p = map(TrialFunction, W)
        v1, v2, q = map(TestFunction, W)

        dxGamma = Measure('dx', domain=iface_mesh)
        # We will need traces of the functions on the boundary
        Tu1, Tu2 = map(lambda x: Trace(x, iface_mesh), (u1, u2))
        Tv1, Tv2 = map(lambda x: Trace(x, iface_mesh), (v1, v2))

        a00 = inner(Constant(sigma['i'])*grad(u1), grad(v1))*dx
        a01 = 0
        a02 = inner(p, Tv1)*dxGamma
        
        a10 = 0
        a11 = inner(Constant(sigma['e'])*grad(u2), grad(v2))*dx 
        a12 = -inner(p, Tv2)*dxGamma

        a20 = inner(q, Tu1)*dxGamma
        a21 = -inner(q, Tu2)*dxGamma
        a22 = -Constant(beta)*inner(p, q)*dxGamma

        a = [[a00, a01, a02], [a10, a11, a12], [a20, a21, a22]]

        # Get the block possibly in terms of products of matrices
        A = ii_assemble(a)
        # As block stuff
        A = ii_convert(A, algorithm='')

        # Boundary conditions for the problem
        Ve_bcs = [DirichletBC(V2, value, exterior_facet_f, tag) for tag, value in bcs.items()]
        bcs = [[], Ve_bcs, []]

        bcs = block_bc(bcs, symmetric=True)
        rhs_bcs = bcs.apply(A)

        # Let's have a function that gives as the rhs of the linear system.
        # Note that forms in L use values of the rhs dictionary so it is
        # through changing these that different b can be obtained. Similarly
        # for the boundary conditions and their bcs dict
        get_L = lambda: [inner(rhs['i'], v1)*dx, inner(rhs['e'], v2)*dx, inner(rhs['m'], q)*dxGamma]

        def get_b(L=get_L, rhs_bcs=rhs_bcs):
            '''Forcing vector with boundary conditions'''
            # Convert to vec from possible mat-vec products
            b = block_vec(map(ii_convert, ii_assemble(L())))
            # Apply bcs to block_vec
            rhs_bcs.apply(b)
            # Make monolithic for the solver
            b = ii_convert(b)

            return b
    
        A = ii_convert(A)  # Monolithic system for the solver
        # Check that this is not a 
        assert not check_eigvals or np.min(np.abs(np.linalg.eigvalsh(A.array()))) > 1E-10

        self.A = A
        self.W = W
        self.__get_b__ = get_b

    @property
    def b(self): return self.__get_b__()

# --------------------------------------------------------------------

if __name__ == '__main__':
    import sympy as sym

    set_log_level(WARNING)
    
    # Sample use test
    if True:
        mesh = UnitSquareMesh(8, 8)

        sigma = {'e': Constant(1), 'i': Constant(1)}
        f = {'e': Constant(1), 'i': Constant(1), 'm': Constant(1)}
        beta = Constant(1E-3)
        bcs = {1: Constant(1), 2: Constant(1), 3: Constant(1), 4: Constant(1)}

        problem = EmiBvpProblem(mesh, sigma, f, beta, bcs, check_eigvals=True)
        A, b, W = problem.A, problem.b, problem.W
        # For solving the block components are converted to a monolithic matrix
        wh = ii_Function(W)
        LUSolver('umfpack').solve(A, wh.vector(), b)
        
        # All the coeffcients should stay away from NaNs and INFs
        is_sane = lambda v: not np.any(np.isnan(v)) and not np.any(np.isinf(v))
        
        assert all(is_sane(x.vector().get_local())  for x in wh)

    # Run MMS for the step solver
    x, y = sym.symbols('x[0], x[1]')
    ui = sym.sin(2*pi*x)*sym.sin(2*pi*y)

    ue = 1 + sym.sin(2*pi*x)*sym.sin(2*pi*y)
    du = ui - ue
    # Laplacians
    fi = -(ui.diff(x, 2) + ui.diff(y, 2))
    fe = -(ue.diff(x, 2) + ue.diff(y, 2))

    ui_exact = Expression(sym.printing.ccode(ui), degree=4, t=0)
    fi = Expression(sym.printing.ccode(fi), degree=4, t=0)

    ue_exact = Expression(sym.printing.ccode(ue), degree=4)
    fe = Expression(sym.printing.ccode(fe), degree=4)

    du_exact = Expression(sym.printing.ccode(du), degree=4, t=0)
    # Constant with declared high degree to shut up errornorm warnings
    p_exact = Expression('0', degree=4)
    
    rhs = {'e': fe, 'i': fi, 'm': du_exact}
    sigma = {'e': Constant(1), 'i': Constant(1)}
    beta = Constant(0)
    bcs = {k: ue_exact for k in (1, 2, 3, 4)}

    msg = '%.2E %g'
    msg = ' '.join([msg] + ['%.2E(%.2f)']*3)

    errors0, h0, table = None, None, []
    for n in (8, 16, 32, 64, 128, 256):
        mesh = UnitSquareMesh(n, n)
        
        problem = EmiBvpProblem(mesh, sigma, rhs, beta, bcs)
        A, b, W = problem.A, problem.b, problem.W
        
        wh = ii_Function(W)
        LUSolver('umfpack').solve(A, wh.vector(), b)
        info('|Ax-b| = %g' % (A*wh.vector() - b).norm('l2'))
        
        uih, ueh, ph = wh
        
        h = ph.function_space().mesh().hmin()

        errors = np.array([errornorm(ui_exact, uih, 'H1'),
                           errornorm(ue_exact, ueh, 'H1'),
                           errornorm(p_exact, ph, 'L2')])
        
        if errors0 is not None:
            rates = np.log(errors/errors0)/np.log(h/h0)
        else:
            rates = [np.nan]*len(errors)
        errors0, h0 = errors, h

        dimW = sum(sub.dim() for sub in wh.function_space())
        row = tuple([h, dimW] + sum(map(list, zip(errors, rates)), []))

        table.append(msg % row)
        print table[-1]
    # Summary
    print '\n'+'\t'.join(['h', 'dofs', '|e_ui|', '|e_ue|', '|e_p|'])
    for row in table: print row
