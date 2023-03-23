from _common import *
from grid import *
from dof_handler import *

print_testing_statement("DoFHandler")
m = Grid(0,1,10)
dh = DoFHandler(m, 2)
print("Number of dofs", dh.n_dofs)