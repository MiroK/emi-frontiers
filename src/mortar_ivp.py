# Mortar FEM for MMS problem from S3.2 of article
#   https://www.frontiersin.org/articles/10.3389/fphy.2017.00048/full
from dolfin import (Expression, Constant, pi, UnitSquareMesh, interpolate,
                    CompiledSubDomain, MeshFunction, FunctionSpace,
                    sqrt, SubMesh, LUSolver)
from xii.assembler.trace_matrix import trace_mat_no_restrict
import dolfin as df
from mortar_bvp import EmiBvpProblem
from scipy.linalg import eigh
from tqdm import tqdm
import sympy as sym
import numpy as np
from xii import *


# Define analytical solution and sources
x, y, t = sym.symbols('x[0], x[1], t')
ui = (1+sym.exp(-t))*sym.sin(2*pi*x)*sym.sin(2*pi*y)

ue = sym.sin(2*pi*x)*sym.sin(2*pi*y)
du = ui - ue
# Laplacians
fi = -(ui.diff(x, 2) + ui.diff(y, 2))
fe = -(ue.diff(x, 2) + ue.diff(y, 2))

ui_exact = Expression(sym.printing.ccode(ui), degree=4, t=0)
fi = Expression(sym.printing.ccode(fi), degree=4, t=0)

ue_exact = Expression(sym.printing.ccode(ue), degree=4)
fe = Expression(sym.printing.ccode(fe), degree=4)

du_exact = Expression(sym.printing.ccode(du), degree=4, t=0)
# Multiplier as grad(ui).n when evaluated on the membrane is 0
p_exact = Expression('0', degree=4)

# NOTE: the potential difference should in principle be computable by interpolation
# from interior/exterior potentials. However, for this to work set_allow_extrapolate
# flag needs to be set on the potentials. Even then the functions are often
# not interpolated correctly. Therefore we perfom this manually
def get_potential_difference(W):
    '''Precompute restriction matrices for potential difference interpolation'''
    V1, V2, Q = W
    
    T1 = df.PETScMatrix(trace_mat_no_restrict(V1, Q))
    T2 = df.PETScMatrix(trace_mat_no_restrict(V2, Q))
    
    def potential_difference(u1, u2, T1=T1, T2=T2, Q=Q):
        '''Trace(u1 - u2) in Q'''
        d = T1*u1.vector()
        d.axpy(-1., T2*u2.vector())
            
        return Function(Q, d)
    return potential_difference


def solve_ivp(ncells, dt):
    '''
    MMS for EMI problem. The membrane ODE v_t = 1/Cm*(p - I_ion) with 
    I_ion = v (passive model) is discretized as v - v_0 = dt/Cm*(p - v0),
    where v_0 is the previous value of the potential difference. This  
    leads to BVP with beta = dt/Cm and f_m = (1 - beta)*v_0
    '''
    # Parameters of bvp
    sigma = {'e': Constant(1), 'i': Constant(1)}
    f = {'e': fe, 'i': fi, 'm': Constant(0)}
    Cm = Constant(1.)
    beta = Constant(dt/Cm)

    bcs = {k: ue_exact for k in (1, 2, 3, 4)}

    mesh = UnitSquareMesh(ncells, ncells)
    # Assemble the matrix problem only once. Changes in bcs and forcing
    # only effect the right hand side vector
    problem = EmiBvpProblem(mesh, sigma, f, beta, bcs)
    A, W = problem.A, problem.W
    # Factorize once
    solver = LUSolver('umfpack')
    solver.set_operator(A)
    
    v0 = ii_Function(W)[2]
    du_exact.t = 0  # Initial conditions
    v0.assign(interpolate(du_exact, v0.function_space()))
    # The force points to v0, so in the by updating v0 (ui_h - ue_h on
    # the interface) we get update f_m as well
    f_m = Constant(1 - beta)*v0
    f['m'] = f_m

    potential_difference = get_potential_difference(W)
    
    # Given du solve for ui, ui, and p
    time, Tstop = 0, 0.1
    du_exact.t = time

    wh = ii_Function(W)
    for i in tqdm(range(int(Tstop/dt))):
        time += dt
        # Update forces and boundary conditions
        fi.t = time
        fe.t = time
        ue_exact.t = time
        # The new rhs for the linear system
        b = problem.b
        # And the corresponding solution
        solver.solve(wh.vector(), b)
        # Update v0 with the new differene of ui_h and ue_h
        v0.assign(potential_difference(wh[0], wh[1]))

    return wh, time


def errornorm(f, fh, kind, **kwargs):
    '''Dolfin errornorm extended to L^inf norm approx'''
    if kind.lower() == 'linf':
        V = FunctionSpace(fh.function_space().mesh(), 'DG', 1)
        e = interpolate(f, V)
        e.vector().axpy(-1, interpolate(fh, V).vector())
        
        return e.vector().norm('linf')

    return df.errornorm(f, fh, kind, **kwargs)


def string_fmt(data, nnorms):
    '''Pretty print row of the table'''
    return ((','.join(['%.2E, %d'] + ['%.2E(%.2f)']*nnorms)) % data).split(',')


def get_errors(ui_h, ue_h, ph):
    '''Solution errors as in Table 3 in the paper'''
    return np.array([
        sqrt(errornorm(ui_exact, ui_h, 'H1')**2 + errornorm(ue_exact, ue_h, 'H1')**2),
        sqrt(errornorm(ui_exact, ui_h, 'L2')**2 + errornorm(ue_exact, ue_h, 'L2')**2),
        errornorm(p_exact, ph, 'L2'),
        max(errornorm(ui_exact, ui_h, 'linf'), errornorm(ue_exact, ue_h, 'Linf')),
        errornorm(du_exact, du_h, 'linf')])

# ---------------------------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import Function, set_log_level, WARNING

    set_log_level(WARNING)

    # Headers for the table
    table = [['h', 'dim(Wh)', '|u-uh|_1', '|u-uh|_0', '|p-ph|_0', '|u-uh|_oo', '|duh-du|_oo']]

    errors0, h0 = None, None
    for n in [8, 16, 32, 64, 128, 256, 512]:
        wh, time = solve_ivp(ncells=n, dt=1E-2/n)

        # We look at error in u_e, u_i
        [setattr(f, 't', time) for f in (ui_exact, ue_exact, du_exact)]

        ui_h, ue_h, ph = wh

        du_h = get_potential_difference(wh.function_space())(ui_h, ue_h)
        
        h = ph.function_space().mesh().hmin()
        # Norms as in table 3 of the paper
        errors = get_errors(ui_h, ue_h, ph)
        
        if h0 is not None:
            rates = np.log(errors/errors0)/np.log(h/h0)
        else:
            rates = [np.nan]*len(errors)
        h0, errors0 = h, errors

        dimW = sum(sub.dim() for sub in wh.function_space())
        data = sum(map(list, zip(errors, rates)), [h, dimW])
        table.append(string_fmt(tuple(data), len(errors)))
    print
    
    col_width = max(len(word) for row in table for word in row)+2
    table = iter(table)
    # Header
    row = next(table)
    print '|'.join(word.ljust(col_width) for word in row)
    # Separator
    print '|'.join(['-'*col_width]*len(row))
    # Actual data
    for row in table:
        print '|'.join(word.ljust(col_width) for word in row)
