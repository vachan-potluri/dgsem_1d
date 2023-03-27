from common import *

class LegendreBasis:
    # A class to store certain quantities related to legendre polynomials
    def __init__(self, N):
        assert N>=1 and N<=5, "Only 1st to 5th order polynomial bases are supported"
        self.N = N
        self.total_variations = np.zeros(N+1)
        for i_poly in range(1,N+1):
            # TV of 0th Legendre polynomial is 0
            poly_coeffs = np.zeros(i_poly+1)
            poly_coeffs[-1] = 1.0
            poly = np.polynomial.Legendre(poly_coeffs)
            poly_der = poly.deriv(1)
            def integ_func(x):
                return abs(np.polynomial.legendre.legval(x, poly_der.coef))
            integ_samples = np.linspace(-1,1,1000)
            integ_vals = integ_func(integ_samples)
            self.total_variations[i_poly] = np.trapz(integ_vals, integ_samples)