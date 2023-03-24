from _common import *
from euler import *

print_testing_statement("Euler")
prim = [1.4,0.3,1]
cons = Euler.prim_to_cons(prim)
print("Conservative state")
print(cons)
prim_aux = Euler.cons_to_prim_aux(cons)
print("Auxiliary primitive state")
print(prim_aux)
print("Exact flux")
print(Euler.get_flux(cons))
print("Chandrashekhar volume flux")
print(Euler.chandrashekhar_volume_flux(cons, cons))
print("K and Kinv")
print(Euler.get_K(prim[1], prim_aux[3], (cons[2]+prim[2])/cons[0]))
print(Euler.get_Kinv(prim[1], prim_aux[3], (cons[2]+prim[2])/cons[0]))
print("Chandrashekhar surface flux")
print(Euler.chandrashekhar_surface_flux(cons, cons, 0.5))