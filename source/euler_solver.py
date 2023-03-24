from common import *

from cell_vector import *
from change_basis_matrix import *
from dof_handler import *
from dof_vector import *
from euler import *
from grid import *

class EulerSolver:
    # The solver class which coordinates everything
    
    def __init__(self):
        # does nothing, just initialises the instance variables to None
        self.mesh = None
        self.dof_handler = None
        self.states = None
        self.alpha = None
    
    def set_mesh(self, mesh):
        # sets the mesh, and variables that are defined on the mesh
        assert type(mesh) == Grid, "Expected Grid type variable"
        print("\nSetting mesh")
        print(f"\tDomain: [{mesh.x_faces[0]}, {mesh.x_faces[-1]}]")
        print(f"\t# cells: {mesh.n_cells}")
        self.mesh = mesh
        print("Initialising blender vector")
        self.alpha = CellVector(mesh)
    
    def distribute_dofs(self, N):
        # modifies the dof handler and solution vectors
        assert N>=1, "Expected N>=1"
        assert self.mesh != None, "Set the mesh before distributing dofs"
        print("\nDistributing DoFs")
        self.dof_handler = DoFHandler(self.mesh, N)
        print(f"\t# dofs: {self.dof_handler.n_dofs}")
        print("Initialising state vector")
        self.states = DoFVector(self.dof_handler, Euler.n_vars)
    
    def set_states(self, func, prim=True):
        # sets the states based on the function provided
        # the function should take position as the input and return the euler state
        # by default, the function is assumed to return primitive variables (rho, u, p)
        pass