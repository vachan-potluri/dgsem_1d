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
            (rhoinv*rhoe*(Euler.gamma-1)*Euler.gamma)**0.5
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
    def get_K(u, a, H):
        # returns the right eigenvector matrix
        K = np.zeros((Euler.n_vars, Euler.n_vars))
        # setting columns
        K[:,0] = [1, u-a, H-u*a]
        K[:,1] = [1, u, 0.5*u**2]
        K[:,2] = [1, u+a, H+u*a]
        return K
    
    @staticmethod
    def get_Kinv(u, a, H):
        # returns the inverse of right eigenvector matrix
        Kinv = np.zeros((Euler.n_vars, Euler.n_vars))
        # setting rows
        temp = 1/(Euler.gamma-1)
        Kinv[0,:] = [H+a*(u-a)*temp, -u-a*temp, 1]
        Kinv[1,:] = [2*(2*a**2*temp - H), 2*u, -2]
        Kinv[2,:] = [H-a*(u+a)*temp, -u+a*temp, 1]
        return Kinv
    
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
        u_avg = 0.5*(priml[1]+primr[1])
        p_hat = 0.5*(consl[0]+consr[0])/(betal+betar)

        flux = np.zeros(Euler.n_vars)
        flux[0] = rho_ln*u_avg
        flux[1] = flux[0]*u_avg + p_hat
        flux[2] = flux[0]*(
            0.5/(beta_ln*(Euler.gamma-1)) + u_avg**2 - 0.25*(priml[1]**2+primr[1]**2)
        ) + p_hat*u_avg
        return flux
    
    @staticmethod
    def chandrashekhar_surface_flux(consl, consr, flux_blender_value=0):
        # chandrashekhar's surface flux
        # doesn't assert positivity

        # first set the volume flux
        priml = Euler.cons_to_prim(consl)
        primr = Euler.cons_to_prim(consr)
        betal = 0.5*consl[0]/priml[2]
        betar = 0.5*consr[0]/primr[2]
        rho_ln = Euler.log_avg(consl[0], consr[0])
        beta_ln = Euler.log_avg(betal, betar)
        u_avg = 0.5*(priml[1]+primr[1])
        p_hat = 0.5*(consl[0]+consr[0])/(betal+betar)

        flux = np.zeros(Euler.n_vars)
        flux[0] = rho_ln*u_avg
        flux[1] = flux[0]*u_avg + p_hat
        flux[2] = flux[0]*(
            0.5/(beta_ln*(Euler.gamma-1)) + u_avg**2 - 0.25*(priml[1]**2+primr[1]**2)
        ) + p_hat*u_avg

        # hybrid matrix based stabilisation
        rho_sql = consl[0]**0.5
        rho_sqr = consr[0]**0.5
        temp = 1/(rho_sql + rho_sqr)
        ui = (rho_sql*priml[1] + rho_sqr*primr[1])*temp # 'i' for intermediate
        Hi = ( (consl[2]+priml[2])/rho_sql + (consr[2]+primr[2])/rho_sqr )*temp
        ai = ((Euler.gamma-1)*(Hi-0.5*ui**2))**0.5
        ui_abs = abs(ui)
        lambda_max_rus = ui_abs + ai
        lambda2_blend = flux_blender_value*lambda_max_rus + (1-flux_blender_value)*ui_abs
        eigenvalues = np.array([lambda_max_rus, lambda2_blend, lambda_max_rus])
        K = Euler.get_K(ui, ai, Hi)
        Kinv = Euler.get_Kinv(ui, ai, Hi)
        A = K.dot( np.diag(eigenvalues).dot(Kinv) )
        flux -= 0.5*A.dot(consr-consl)
        return flux
    