import sys
sys.path.append("../../../source/")

from euler_solver import *

s = EulerSolver()
m = Grid(0,1,100)
s.set_mesh(m)
s.distribute_dofs(3)
t_end = 0.035 # end time
def ic_func_prim(x_dof, x_cell):
    if x_cell<0.3: return np.array([5.99924,19.5975,460.894])
    else: return np.array([5.99242,-6.19633,46.095])
s.set_states(ic_func_prim)
s.set_surface_flux()
s.set_volume_flux()
s.set_blender_params()
def bc_left(t,cons):
    return Euler.prim_to_cons(np.array([5.99924,19.5975,460.894]))
def bc_right(t,cons):
    return Euler.prim_to_cons(np.array([5.99242,-6.19633,46.095]))
s.set_bc_funcs(bc_left, bc_right)
s.set_time_controls(0,t_end,CFL=0.5,rk_order=3)
s.set_write_controls(50)
s.run()