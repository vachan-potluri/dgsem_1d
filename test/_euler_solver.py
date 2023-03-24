from _common import *
from euler_solver import *

s = EulerSolver()
m = Grid(0,1,10)
s.set_mesh(m)
s.distribute_dofs(3)