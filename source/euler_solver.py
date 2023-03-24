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
        self.alpha = None
        self.bc_func_left = None
        self.bc_func_right = None
        self.cbm = None
        self.dof_handler = None
        self.mesh = None
        self.states = None
        self.surface_flux = None
        self.volume_flux = None
    
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
        print("Setting change of basis matrix")
        self.cbm = ChangeBasisMatrix(N)
    
    def set_states(self, func, prim=True):
        # sets the states based on the function provided
        # the function should take dof location and cell center location as the input,
        # and return the euler state
        # by default, the function is assumed to return primitive variables (rho, u, p)
        assert self.states != None, "Dofs must be distributed before setting states"
        for i_cell in range(self.mesh.n_cells):
            for i_dof in range(i_cell*(self.dof_handler.N+1), (i_cell+1)*(self.dof_handler.N+1)):
                if prim==True:
                    self.states.entries[i_dof] = Euler.prim_to_cons(
                        func(
                            self.dof_handler.x_dofs[i_dof],
                            self.mesh.x_cells[i_cell]
                        )
                    )
                else:
                    self.states.entries[i_dof] = func(
                        self.dof_handler.x_dofs[i_dof],
                        self.mesh.x_cells[i_cell]
                    )
    
    def set_bc_funcs(self, left_func, right_func):
        # sets the BC functions
        # the functions should take location, time and inner state as input,
        # and return the ghost state
        self.bc_func_left = left_func
        self.bc_func_right = right_func
    
    def set_surface_flux(self, name="chandrashekhar"):
        # sets the surface flux function
        if name == "chandrashekhar":
            self.surface_flux = Euler.chandrashekhar_surface_flux
        else:
            assert False, "Currently only chandrashekhar flux supported"
    
    def set_volume_flux(self, name="chandrashekhar"):
        # sets the volume flux function
        if name == "chandrashekhar":
            self.volume_flux = Euler.chandrashekhar_volume_flux
        else:
            assert False, "Currently only chandrashekhar flux supported"
    
    def calc_blender(self):
        # calculates the blender values
        # currently uses "pressure times density"
        # currently uses Hennemenn's algorithm (as opposed to Persson's)
        Np1 = self.dof_handler.N+1
        pxrho_nodal_values = np.zeros(Np1)
        threshold = 0.5*10**(-1.8*Np1**0.25)

        # first loop: calculate and clip
        for i_cell in range(self.mesh.n_cells):
            for i_dof_local in range(Np1):
                prim = Euler.cons_to_prim(self.states.entries[i_cell*(Np1)+i_dof_local])
                pxrho_nodal_values[i_dof_local] = prim[0]*prim[2]
            pxrho_modes = self.cbm.get_modes(pxrho_nodal_values)
            energies = pxrho_modes*pxrho_modes
            trouble = max(energies[-1]/np.sum(energies), energies[-2]/np.sum(energies[:-1]))
            alpha = 1/(1+np.exp(-9.21024*(trouble/threshold-1)))
            if alpha < 1e-3:
                self.alpha.entries[i_cell] = 0.0
            elif alpha > 0.5:
                self.alpha.entries[i_cell] = 0.5
            else:
                self.alpha.entries[i_cell] = alpha
        
        # second loop: diffuse
        alpha_old = self.alpha.entries.copy()
        for i_cell in range(self.mesh.n_cells):
            if i_cell == 0:
                i_nbs = [i_cell+1]
            elif i_cell == self.mesh.n_cells-1:
                i_nbs = [i_cell-1]
            else:
                i_nbs = [i_cell-1, i_cell+1]
            self.alpha.entries[i_cell] = max(
                alpha_old[i_cell],
                0.5*np.max(alpha_old[i_nbs])
            )
    
    def calc_rhs(self):
        # calculates the rhs
        rhs = np.zeros(self.dof_handler.n_dofs, dtype=object)
        N = self.dof_handler.N
        
        # first assert positivity
        for i_dof in range(self.dof_handler.n_dofs):
            Euler.assert_positivity(self.states.entries[i_dof])
        
        # calculate blender
        self.calc_blender()

        # calculate the surface (cell interface/boundary) fluxes
        surface_fluxes = np.zeros(self.mesh.n_cells+1, dtype=object)
        cons = self.states.entries[0]
        surface_fluxes[0] = self.surface_flux(
            self.bc_func_left(self.mesh.x_faces[0], 0, cons), # <-- time controller?
            cons,
            self.alpha.entries[0]
        )
        cons = self.states.entries[-1]
        surface_fluxes[-1] = self.surface_flux(
            cons,
            self.bc_func_right(self.mesh.x_faces[-1], 0, cons),
            self.alpha.entries[-1]
        )
        for i_face in range(1, self.mesh.n_cells):
            i_dof_left = (i_face-1)*(N+1) + N
            i_dof_right = i_face*(N+1)
            surface_fluxes[i_face] = self.surface_flux(
                self.states.entries[i_dof_left],
                self.states.entries[i_dof_right],
                0.5*(self.alpha.entries[i_face-1] + self.alpha.entries[i_face])
            )
        
        # calculate the DG residual
        for i_cell in range(self.mesh.n_cells):
            local_dofs = np.arange(i_cell*(N+1), (i_cell+1)*(N+1))
            rhs[local_dofs] = np.zeros(Euler.n_vars)
            # calculate volume flux at every flux point
            