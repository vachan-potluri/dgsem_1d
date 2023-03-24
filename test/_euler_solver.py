from _common import *
from euler_solver import *

s = EulerSolver()
m = Grid(-1,2,10)
s.set_mesh(m)
s.distribute_dofs(3)
def ic_func_prim(x_dof, x_cell):
    if x_cell<0: return np.array([1,0.75,1])
    else: return np.array([0.125,0,0.1])
s.set_states(ic_func_prim)
for i in range(s.dof_handler.n_dofs):
    print(f"x = {s.dof_handler.x_dofs[i]:1.4f}, state = {s.states.entries[i]}")