import sys
sys.path.append("../../source/")

from euler_solver import *

s = EulerSolver()
m = Grid(0,1,67)
s.set_mesh(m)
s.distribute_dofs(2)
def ic_func_prim(x_dof, x_cell):
    if x_cell<0.3: return np.array([1,0.75,1])
    else: return np.array([0.125,0,0.1])
s.set_states(ic_func_prim)
s.set_surface_flux()
s.set_volume_flux()
def bc_left(t,cons):
    return Euler.prim_to_cons(np.array([1,0.75,1]))
def bc_right(t,cons):
    return Euler.prim_to_cons(np.array([0.125,0,0.1]))
s.set_bc_funcs(bc_left, bc_right)
s.set_time_controls(0,0.2,10,CFL=0.1)
s.run()