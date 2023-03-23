from common import *
from q_gauss_lobatto import *
class ChangeBasisMatrix:
    # A class to calculate and store a matrix corresponding to basis conversion from nodal
    # interpolation in Gauss-Lobatto points to modal representation using Legendre series
    # The matrix multiplied by the nodal values will give the modes
    def __init__(self, N):
        # initialises the class with polynomial degree
        assert N>=1 and N<=5, "Only 1st to 5th polynomial basis degrees are supported"
        self.matrix = np.zeros((N+1, N+1))
        gl_points = QGaussLobatto(N+1).q_points
        for col in range(N+1):
            # roots for col-th nodal basis function
            nodal_fn_roots = np.delete(gl_points, col)
            # coefficients, scaled below
            nodal_fn_coeffs = np.polynomial.polynomial.polyfromroots(nodal_fn_roots)
            for root in nodal_fn_roots:
                nodal_fn_coeffs /= (gl_points[col] - root)
            
            # convert to modal basis
            self.matrix[:, col] = np.polynomial.legendre.poly2leg(nodal_fn_coeffs)
    
    def get_modes(self, nodal_values):
        return self.matrix.dot(nodal_values)
