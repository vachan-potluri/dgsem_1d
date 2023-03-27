from common import *
from dof_handler import *

class DoFVector:
    # Class for storing and operating on a solution vector defined on dofs
    # the storage is like an array of structs
    def __init__(self, dh, n_vars=1):
        assert type(dh) == DoFHandler, "Expected a DoFHandler instance as input"
        assert n_vars>=1, "Zero/negative number of variables received"
        self.dof_handler = dh
        if n_vars == 1:
            self.entries = np.zeros(dh.n_cells, dtype=float)
        else:
            self.entries = np.zeros(dh.n_dofs, dtype=object)
            for i in range(dh.n_dofs):
                self.entries[i] = np.zeros(n_vars)
    
    def sample_in_cell(self, i_cell, n_samples):
        # returns the solution evaluated at `n_samples` uniformly spaced sampling points in cell
        # with index `i_cell`
        # the left and right boundaries are included
        # the interpolation is done in reference space
        sampled_solns = np.zeros(n_samples, dtype=object)
        for i in range(n_samples):
            sampled_solns[i] = np.zeros_like(self.entries[i])
        xi_samples = np.linspace(-1,1,n_samples) # sampling points in reference space
        N = self.dof_handler.N
        for i_poly in range(N+1):
            poly_roots = np.delete(self.dof_handler.quad.q_points, i_poly)
            poly_coeffs = np.polynomial.polynomial.polyfromroots(poly_roots)
            for root in poly_roots:
                poly_coeffs /= (self.dof_handler.quad.q_points[i_poly] - root)
            for i_sample, xi_sample in enumerate(xi_samples):
                sampled_solns[i_sample] += (
                    self.entries[i_cell*(N+1)+i_poly]*np.polynomial.polynomial.polyval(
                        xi_sample,
                        poly_coeffs
                    )
                )
        return sampled_solns