from common import *

class Euler:
    # A class for a collection of functions
    # Provides some useful functions to deal with euler states
    gamma = 1.4
    R = 286.7
    n_vars = 3

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
        prim = Euler.cons_to_prim(cons)
        return np.array([cons[1], cons[1]*prim[1]+prim[2], (cons[2]+prim[2])*prim[1]])

    @staticmethod
    def prim_to_cons(prim):
        # doesn't assert positivity
        return np.array([prim[0], prim[0]*prim[1], prim[2]/(Euler.gamma-1)+prim[0]*(prim[1]**2)])
    
    @staticmethod
    def log_avg(x1, x2):
        # to calculate logarithmic average in a numerically stable way
        xi = x1/x2
        f = (xi-1)/(xi+1)
        u = f*f
        if u<1e-2:
            F = 1 + u/3 + u**2/5 + u**3/7
        else:
            F = 0.5*np.log(xi)/f
        return 0.5*(x1 + x2)/F
    
    @staticmethod
    def chandrashekhar_volume_flux(consl, consr):
        # chandrashekhar's volume flux
        # doesn't assert positivity
        priml = Euler.cons_to_prim(consl)
        primr = Euler.cons_to_prim(consr)
        betal = 0.5*consl[0]/priml[2]
        betar = 0.5*consr[0]/primr[2]
        rho_ln = Euler.log_avg(consl[0], consr[0])
        beta_ln = Euler.log_avg(betal, betar)
        rho_avg = 0.5*(consl[0]+consr[0])
        beta_avg = 0.5*(betal+betar)
        u_avg = 0.5*(priml[1]+primr[1])
        p_hat = 0.5*rho_avg/beta_avg

        flux = np.zeros(Euler.n_vars)
        flux[0] = rho_ln*u_avg
        flux[1] = flux[0]*u_avg + p_hat
        flux[2] = flux[0]*(
            0.5/(beta_ln*(Euler.gamma-1)) + u_avg**2 - 0.5*(priml[1]**2+primr[1]**2)
        ) + p_hat*u_avg
        return flux
    