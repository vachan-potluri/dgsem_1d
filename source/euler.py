from common import *

class Euler:
    # A class for a collection of functions
    # Provides some useful functions to deal with euler states
    gamma = 1.4
    R = 286.7

    @staticmethod
    def assert_positivity(cons):
        # asserts positivity
        assert cons[0]>0, "Negative density"
        rhoinv = 1/cons[0]
        rhoe = cons[2] - 0.5*(cons[1]**2)*rhoinv
        assert rhoe>0, "Negative energy"

    @staticmethod
    def cons_to_prim(cons):
        # converts conservative state to primitive state (rho, u, p)
        # doesn't assert positivity
        rhoinv = 1/cons[0]
        rhoe = cons[2] - 0.5*(cons[1]**2)*rhoinv
        return np.array([cons[0], cons[1]*rhoinv, (Euler.gamma-1)*rhoe])
    
    @staticmethod
    def cons_to_prim_aux(cons):
        # converts to primitive variables, and additionally also returns the speed of sound
        # doesn't assert positivity
        rhoinv = 1/cons[0]
        rhoe = cons[2] - 0.5*(cons[1]**2)*rhoinv
        return np.array([
            cons[0],
            cons[1]*rhoinv,
            (Euler.gamma-1)*rhoe,
            (rhoe*(Euler.gamma-1)*Euler.gamma)**0.5
        ])

    @staticmethod
    def get_flux(cons):
        # returns the analytical flux
        prim = cons_to_prim(cons)
        return np.array([cons[1], cons[1]*prim[1]+prim[2], (cons[2]+prim[2])*prim[1]])

    @staticmethod
    def prim_to_cons(prim):
        # doesn't assert positivity
        return np.array([prim[0], prim[0]*prim[1], prim[2]/(Euler.gamma-1)+prim[0]*(prim[1]**2)])

        
    
    
    