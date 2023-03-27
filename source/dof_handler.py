from common import *
from grid import *
from q_gauss_lobatto import *

class DoFHandler:
    # A class to distribute DoFs on a 1d mesh
    def __init__(self, mesh, N):
        assert N>=1 and N<=5, "Currently only values of N supported are [1,5]"
        assert type(mesh) == Grid, "Invalid mesh input. Expected a Grid class instance"
        self.N = N
        self.mesh = mesh
        self.n_dofs = mesh.n_cells*(N+1)
        self.quad = QGaussLobatto(N+1)
        self.x_dofs = np.zeros(self.n_dofs)
        for i_cell in range(mesh.n_cells):
            self.x_dofs[i_cell*(N+1):(i_cell+1)*(N+1)] = (
                (self.quad.q_points+1)*mesh.J + mesh.x_faces[i_cell]
            )
        self.D = np.zeros((N+1,N+1)) # polynomial derivative matrix
        for col in range(N+1):
            # construct col-th basis polynomial in reference space
            # and also compute its derivative
            roots = np.delete(self.quad.q_points, col)
            basis_mono_coeffs = np.polynomial.polynomial.polyfromroots(roots)
            for root in roots:
                basis_mono_coeffs /= (self.quad.q_points[col]-root)
            basis_der_mono_coeffs = np.polynomial.polynomial.polyder(basis_mono_coeffs)
            for row in range(N+1):
                self.D[row,col] = np.polynomial.polynomial.polyval(
                    self.quad.q_points[row],
                    basis_der_mono_coeffs
                )