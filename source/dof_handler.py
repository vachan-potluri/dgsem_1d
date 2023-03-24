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