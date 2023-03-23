from _common import *
from euler import *

print_testing_statement("Euler")
prim = [1,0,1]
cons = Euler.prim_to_cons(prim) 
print(cons)
prim_aux = Euler.cons_to_prim_aux(cons)
print(prim_aux)