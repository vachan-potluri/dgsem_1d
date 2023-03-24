from _common import *
from grid import *
from cell_vector import *

print_testing_statement("CellVector")
m = Grid(0,1,10)
soln = CellVector(m, 3)
print(soln.entries)