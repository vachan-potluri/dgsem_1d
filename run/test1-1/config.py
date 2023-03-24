import sys
sys.path.append("../../source/")

from euler_solver import *

s = EulerSolver()
m = Grid(0,1,50)
s.set_mesh(m)
s.distribute_dofs(3)
def ic_func_prim(x_dof, x_cell):
    if x_cell<0.3: return np.array([1,0.75,1])
    else: return np.array([0.125,0,0.1])
s.set_states(ic_func_prim)
s.set_surface_flux()
s.set_volume_flux()
def bc_free(x,t,cons):
    return cons
s.set_bc_funcs(bc_free, bc_free)
s.set_time_controls(0,0.2,1)
s.run()