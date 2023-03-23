from common import *
from grid import *

class DoFHandler:
    # A class to distribute DoFs on a 1d mesh
    def __init__(self, mesh, N):
        assert N>=1 and N<=5, "Currently only values of N supported are [1,5]"
        assert type(mesh) is Grid, "Invalid mesh input. Expected a Grid class instance"
        self.N = N
        self.mesh = mesh
        self.n_dofs = mesh.n_cells*(N+1)