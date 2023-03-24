from common import *

class Grid:
    # A simple class to store and use 1d mesh
    # Currently works for uniform mesh only
    def __init__(self, x_left, length, n_cells):
        assert length>0, "Negative/zero length received"
        assert n_cells>0, "Negative/zero number of cells received"
        self.n_cells = n_cells
        self.x_faces = np.linspace(x_left, x_left+length, n_cells+1) # cell boundaries
        self.dx = self.x_faces[1] - self.x_faces[0] # cell width (uniform)
        self.x_cells = np.linspace(x_left+0.5*self.dx, x_left+length-0.5*self.dx, n_cells)
        self.J = 0.5*self.dx # Jacobian of mapping from any cell to reference space [-1,1]
        self.Jinv = 1.0/self.J