from _common import *
from grid import *
from dof_handler import *

print_testing_statement("DoFHandler")
m = Grid(-1,4,2)
dh = DoFHandler(m, 3)
print("Number of dofs", dh.n_dofs)
print("DoF locations")
print(dh.x_dofs)