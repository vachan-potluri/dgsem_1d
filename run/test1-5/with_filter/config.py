import sys
sys.path.append("../../../source/")

from euler_solver import *

s = EulerSolver()
m = Grid(0,1,100)
s.set_mesh(m)
s.distribute_dofs(3)
def ic_func_prim(x_dof, x_cell):
    if x_cell<0.3: return np.array([1,-19.59745,1000])
    else: return np.array([1,-19.59745,0.01])
s.set_states(ic_func_prim)
s.set_surface_flux()
s.set_volume_flux()
s.set_blender_params(do_filter_noise=True)
def bc_left(t,cons):
    return Euler.prim_to_cons(np.array([1,-19.59745,1000]))
def bc_right(t,cons):
    return Euler.prim_to_cons(np.array([1,-19.59745,0.01]))
s.set_bc_funcs(bc_left, bc_right)
s.set_time_controls(0,0.012,CFL=0.5,rk_order=3)
s.set_write_controls(50)
s.run()