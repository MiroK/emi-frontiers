from dolfin import (Expression, Constant, pi, UnitSquareMesh, interpolate,
                    errornorm, CompiledSubDomain, MeshFunction, FunctionSpace,
                    sqrt, SubMesh, LUSolver)
from mortar_bvp import EmiBvpProblem
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


def solve_ivp(ncells, dt):
    '''MMS for EMI problem'''
    # Parameters of bvp
    sigma = {'e': Constant(1), 'i': Constant(1)}
    f = {'e': fe, 'i': fi, 'm': Constant(0)}
    Cm = Constant(1.)
    beta = Constant(Cm/dt)  # Cm/dt
    bcs = {k: ue_exact for k in (1, 2, 3, 4)}

    mesh = UnitSquareMesh(ncells, ncells)
    problem = EmiBvpProblem(mesh, sigma, f, beta, bcs)
    A, W = problem.A, problem.W
    
    solver = LUSolver('umfpack')
    solver.set_operator(A)
    
    wh = ii_Function(W)
    p0 = wh[2]
    du_exact.t = 0
    p0.assign(interpolate(du_exact, p0.function_space()))
    # Update dict
    f_m = (1-beta**(-1))*p0
    f['m'] = f_m
    

    # Given du solve for ui, ui, and p
    time, Tstop = 0, 0.1
    du_exact.t = time

    for i in tqdm(range(int(Tstop/dt))):
        time += dt
        # Update forces and boundary conditions
        fi.t = time
        fe.t = time
        ue_exact.t = time
        # The new rhs for the linear system
        b = problem.b
        # NOTE: 
        solver.solve(wh.vector(), b)
    
    time = time - dt

    return wh, time

# ---------------------------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import File, set_log_level, WARNING

    set_log_level(WARNING)

    errors0, h0, table = None, None, []
    for n in [16, 32, 64, 128, 256]:
        print n
        wh, time = solve_ivp(ncells=n, dt=0.01/n)

        # We look at error in u_e, u_i
        [setattr(f, 't', time) for f in (ui_exact, ue_exact, du_exact)]

        ui_h, ue_h, ph = wh
        h = ph.function_space().mesh().hmin()
        errors = np.array([errornorm(ui_exact, ui_h, 'H1'),
                           errornorm(ue_exact, ue_h, 'H1'),
                           errornorm(p_exact, ph, 'L2')])

        if h0 is not None:
            rates = np.log(errors/errors0)/np.log(h/h0)
        else:
            rates = [np.nan]*len(errors)
        h0, errors0 = h, errors
        
        data = sum(map(list, zip(errors, rates)), [n, h])

        dimW = sum(sub.dim() for sub in wh.function_space())
        row = tuple([h, dimW] + sum(map(list, zip(errors, rates)), []))

        msg = ' | '.join(['%.2E %d'] + ['%.3E(%.2f)']*len(errors))
        table.append(msg % row)
        print table[-1]

    # Summary
    print '\n'+'\t'.join(['h', 'dofs', '|e_ui|', 'r', '|e_ue|', 'r', '|e_p|', 'r'])
    for row in table: print row

    File('foo.pvd') << ph
