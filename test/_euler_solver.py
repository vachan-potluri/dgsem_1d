from _common import *
from euler_solver import *

s = EulerSolver()
m = Grid(-1,2,10)
s.set_mesh(m)
s.distribute_dofs(3)
def ic_func_prim(x_dof, x_cell):
    if x_dof<0.1: return np.array([1,0.75,1])
    else: return np.array([0.125,0,0.1])
s.set_states(ic_func_prim)
for i in range(s.dof_handler.n_dofs):
    print(f"x = {s.dof_handler.x_dofs[i]:1.4f}, state = {s.states.entries[i]}")
s.set_surface_flux()
s.set_volume_flux()
def bc_free(x,t,cons):
    return cons
s.set_bc_funcs(bc_free, bc_free)

s.calc_blender()
print("Blender values:")
print(s.alpha.entries)
s.calc_rhs()