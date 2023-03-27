from common import *
from grid import *

class CellVector:
    # Class to store a solution vector that is defined on cells rather than dofs
    # the storage is like an array of structs
    def __init__(self, mesh, n_vars=1):
        assert type(mesh) == Grid, "Expected a DoFHandler instance as input"
        assert n_vars>=1, "Zero/negative number of variables received"
        self.mesh = mesh
        if n_vars == 1:
            self.entries = np.zeros(mesh.n_cells, dtype=float)
        else:
            self.entries = np.zeros(mesh.n_cells, dtype=object)
            for i in range(mesh.n_cells):
                self.entries[i] = np.zeros(n_vars)