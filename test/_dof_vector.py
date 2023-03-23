from _common import *
from grid import *
from dof_handler import *
from dof_vector import *

print_testing_statement("DoFVector")
m = Grid(0,1,10)
dh = DoFHandler(m, 2)
soln = DoFVector(dh, 3)
print(soln.entries)