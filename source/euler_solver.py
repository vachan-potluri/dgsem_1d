from common import *

from cell_vector import *
from change_basis_matrix import *
from dof_handler import *
from dof_vector import *
from euler import *
from grid import *
from legendre_basis import *
import ls3_rk45_coeffs

class EulerSolver:
    # The solver class which coordinates everything
    
    def __init__(self):
        # does nothing, just initialises the instance variables to None
        self.alpha = None
        self.alpha_max = None
        self.do_calc_noise = None
        self.do_filter_noise = None
        self.bc_func_left = None
        self.bc_func_right = None
        self.cbm = None
        self.CFL = None
        self.D = None
        self.dof_handler = None
        self.do_time_step = None
        self.legendre_basis = None
        self.mesh = None
        self.n_samples_per_cell = None
        self.states = None
        self.surface_flux = None
        self.time = None
        self.end_time = None
        self.volume_flux = None
        self.write_freq = None
    
    def set_mesh(self, mesh):
        # sets the mesh, and variables that are defined on the mesh
        assert type(mesh) == Grid, "Expected Grid type variable"
        print("\nSetting mesh")
        print(f"\tDomain: [{mesh.x_faces[0]}, {mesh.x_faces[-1]}]")
        print(f"\t# cells: {mesh.n_cells}")
        self.mesh = mesh
        print("Initialising blender vector")
        self.alpha = CellVector(mesh, 2) # first variable: shock alpha, second variable: noise indicator
    
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
        print("Setting legendre basis")
        self.legendre_basis = LegendreBasis(N)
    
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
        # the functions should take time and inner state as input,
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
    
    def set_blender_params(
        self, alpha_min=1e-3, alpha_max=0.5, do_calc_noise=True, do_filter_noise=False
    ):
        self.alpha_min = alpha_min
        self.alpha_max = alpha_max
        self.do_calc_noise = do_calc_noise
        self.do_filter_noise = False
        if do_calc_noise:
            self.do_filter_noise = do_filter_noise
    
    def set_time_controls(self, start_time, end_time, CFL=0.5, rk_order=3):
        self.time = start_time
        self.end_time = end_time
        self.CFL = CFL
        if rk_order == 3:
            self.do_time_step = self.tvd_rk3_time_step
        elif rk_order == 4:
            self.do_time_step = self.ls3_rk45_time_step
        else:
            assert False, "Currently only 3rd and 4th order RK integration are supported"
    
    def set_write_controls(self, write_freq=50, n_samples_per_cell=15):
        self.write_freq = write_freq
        self.n_samples_per_cell = n_samples_per_cell
    
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
            if alpha < self.alpha_min:
                self.alpha.entries[i_cell][0] = 0.0
            elif alpha > self.alpha_max:
                self.alpha.entries[i_cell][0] = self.alpha_max
            else:
                self.alpha.entries[i_cell][0] = alpha
        
        # second loop: diffuse
        alpha_old = np.zeros(self.mesh.n_cells)
        for i_cell in range(self.mesh.n_cells):
            alpha_old[i_cell] = self.alpha.entries[i_cell][0]
        for i_cell in range(self.mesh.n_cells):
            if i_cell == 0:
                i_nbs = [i_cell+1]
            elif i_cell == self.mesh.n_cells-1:
                i_nbs = [i_cell-1]
            else:
                i_nbs = [i_cell-1, i_cell+1]
            self.alpha.entries[i_cell][0] = max(
                alpha_old[i_cell],
                0.5*np.max(alpha_old[i_nbs])
            )
    
    def calc_noise(self):
        # calculates total variation based noise values
        # currently uses "pressure times density"
        Np1 = self.dof_handler.N+1
        pxrho_nodal_values = np.zeros(Np1)
        for i_cell in range(self.mesh.n_cells):
            if self.alpha.entries[i_cell][0] < self.alpha_min:
                for i_dof_local in range(Np1):
                    prim = Euler.cons_to_prim(self.states.entries[i_cell*(Np1)+i_dof_local])
                    pxrho_nodal_values[i_dof_local] = prim[0]*prim[2]
                pxrho_modes = self.cbm.get_modes(pxrho_nodal_values)
                pxrho_abs_modes = abs(pxrho_modes)
                scaled_tvs = pxrho_abs_modes*self.legendre_basis.total_variations
                tv_bound = np.sum(scaled_tvs)
                if tv_bound <= 1e-2*np.sum(pxrho_abs_modes):
                    # negligible noise
                    self.alpha.entries[i_cell][1] = 0.0
                else:
                    self.alpha.entries[i_cell][1] = (
                        (0.5*self.alpha_max - self.alpha.entries[i_cell][0])*
                        (1 - scaled_tvs[1]/tv_bound)
                    )
                if self.do_filter_noise:
                    self.alpha.entries[i_cell][0] += self.alpha.entries[i_cell][1]
            else:
                self.alpha.entries[i_cell][1] = 0.0
    
    def calc_rhs(self):
        # calculates the rhs
        rhs = np.zeros(self.dof_handler.n_dofs, dtype=object)
        N = self.dof_handler.N
        
        # first assert positivity
        for i_dof in range(self.dof_handler.n_dofs):
            Euler.assert_positivity(self.states.entries[i_dof])
        
        # calculate blender
        self.calc_blender()
        if self.do_calc_noise:
            self.calc_noise()

        # calculate the surface (cell interface/boundary) fluxes
        surface_fluxes = np.zeros(self.mesh.n_cells+1, dtype=object)
        cons = self.states.entries[0]
        surface_fluxes[0] = self.surface_flux(
            self.bc_func_left(self.time, cons),
            cons,
            self.alpha.entries[0][0]
        )
        cons = self.states.entries[-1]
        surface_fluxes[-1] = self.surface_flux(
            cons,
            self.bc_func_right(self.time, cons),
            self.alpha.entries[-1][0]
        )
        for i_face in range(1, self.mesh.n_cells):
            i_dof_left = (i_face-1)*(N+1) + N
            surface_fluxes[i_face] = self.surface_flux(
                self.states.entries[i_dof_left],
                self.states.entries[i_dof_left+1],
                0.5*(self.alpha.entries[i_face-1][0] + self.alpha.entries[i_face][0])
            )
        
        # calculate the DG residual
        for i_cell in range(self.mesh.n_cells):
            cell_dofs = np.arange(i_cell*(N+1), (i_cell+1)*(N+1))
            for i_dof in cell_dofs:
                rhs[i_dof] = np.zeros(Euler.n_vars)
            # calculate volume flux at every internal flux point and add its contrib to rhs
            # internal flux points contribution
            for i in range(1,N+1):
                vol_flux = np.zeros(Euler.n_vars)
                for j in range(i,N+1):
                    for k in range(i):
                        vol_flux += (
                            2*self.dof_handler.quad.q_weights[k]*self.dof_handler.D[k,j]*
                            self.volume_flux(
                                self.states.entries[cell_dofs[k]],
                                self.states.entries[cell_dofs[j]]
                            )
                        )
                rhs[cell_dofs[i-1]] += vol_flux/self.dof_handler.quad.q_weights[i-1]
                rhs[cell_dofs[i]] += -vol_flux/self.dof_handler.quad.q_weights[i]
            # surface contribution
            rhs[cell_dofs[0]] += -surface_fluxes[i_cell]/self.dof_handler.quad.q_weights[0]
            rhs[cell_dofs[N]] += surface_fluxes[i_cell+1]/self.dof_handler.quad.q_weights[N]
        
        # update rhs for troubled cells
        for i_cell in range(self.mesh.n_cells):
            if self.alpha.entries[i_cell][0] > 0:
                # scale the rhs
                cell_alpha = self.alpha.entries[i_cell][0]
                cell_dofs = np.arange(i_cell*(N+1), (i_cell+1)*(N+1))
                for i_dof in cell_dofs:
                    rhs[i_dof] *= (1-cell_alpha)
                # internal subcell fluxes' contribution
                for i in range(1,N+1):
                    subcell_surf_flux = self.surface_flux(
                        self.states.entries[cell_dofs[i-1]],
                        self.states.entries[cell_dofs[i]],
                        cell_alpha
                    )
                    rhs[cell_dofs[i-1]] += (
                        cell_alpha*subcell_surf_flux/self.dof_handler.quad.q_weights[i-1]
                    )
                    rhs[cell_dofs[i]] += (
                        -cell_alpha*subcell_surf_flux/self.dof_handler.quad.q_weights[i]
                    )
                # surface contribution
                rhs[cell_dofs[0]] += (
                    -cell_alpha*surface_fluxes[i_cell]/self.dof_handler.quad.q_weights[0]
                )
                rhs[cell_dofs[N]] += (
                    cell_alpha*surface_fluxes[i_cell+1]/self.dof_handler.quad.q_weights[N]
                )
        rhs *= -self.mesh.Jinv
        return rhs
    
    def calc_time_step(self):
        # calculate and return the time step
        N = self.dof_handler.N
        dt = 1e10
        for i_cell in range(self.mesh.n_cells):
            dt_cell = 1e10
            cell_dofs = np.arange(i_cell*(N+1), (i_cell+1)*(N+1))
            for i in range(N+1):
                cons = self.states.entries[cell_dofs[i]]
                prim_aux = Euler.cons_to_prim_aux(cons)
                dt_cell = min(
                    dt_cell,
                    self.CFL/N**1.5*(
                        1/(self.mesh.Jinv*abs(prim_aux[1]) + prim_aux[3]/self.mesh.dx)
                    )
                )
            dt = min(dt, dt_cell)
        return dt
    
    def tvd_rk3_time_step(self, dt):
        # performs TVD-RK3 update with the given time step
        states_old = self.states.entries.copy()
        rhs = self.calc_rhs()
        self.states.entries += dt*rhs
        rhs = self.calc_rhs()
        self.states.entries = 0.75*states_old + 0.25*(self.states.entries + dt*rhs)
        rhs = self.calc_rhs()
        self.states.entries = (states_old + 2*(self.states.entries + dt*rhs))/3.0
        self.time += dt
    
    def ls3_rk45_time_step(self, dt):
        # performs 3-register 5-stage RK4 update
        # stage 1
        rhs = self.calc_rhs()
        self.states.entries += ls3_rk45_coeffs.a_outer[0]*dt*rhs
        rhs1 = rhs.copy()
        # stage 2
        rhs = self.calc_rhs()
        self.states.entries += dt*(
            (ls3_rk45_coeffs.a_inner[0]-ls3_rk45_coeffs.a_outer[0])*rhs1 +
            ls3_rk45_coeffs.a_outer[1]*rhs
        )
        rhs2 = rhs1.copy()
        rhs1 = rhs.copy()
        # stage 3
        for k in [0,1,2]:
            rhs = self.calc_rhs()
            self.states.entries += dt*(
                (ls3_rk45_coeffs.b[k] - ls3_rk45_coeffs.a_inner[k])*rhs2 +
                (ls3_rk45_coeffs.a_inner[k+1] - ls3_rk45_coeffs.a_outer[k+1])*rhs1 +
                ls3_rk45_coeffs.a_outer[k]*rhs
            )
            rhs2 = rhs1.copy()
            rhs1 = rhs.copy()
        self.time += dt
    
    def write_solution(self, counter):
        print("Writing solution")
        filename = f"solution_{counter:06d}.csv"
        save_array = np.zeros((self.mesh.n_cells*self.n_samples_per_cell, 6))
        for i_cell in range(self.mesh.n_cells):
            sampled_states = self.states.sample_in_cell(i_cell, self.n_samples_per_cell)
            x_samples = np.linspace(
                self.mesh.x_faces[i_cell],
                self.mesh.x_faces[i_cell+1],
                self.n_samples_per_cell
            )
            for i_sample, x_sample in enumerate(x_samples):
                prim = Euler.cons_to_prim(sampled_states[i_sample])
                save_array[i_cell*self.n_samples_per_cell + i_sample, :] = [
                    x_sample,
                    prim[0],
                    prim[1],
                    prim[2],
                    self.alpha.entries[i_cell][0],
                    self.alpha.entries[i_cell][1]
                ]
        np.savetxt(filename, save_array, delimiter=",", header="x,rho,u,p,alpha_shock,alpha_noise")
    
    def run(self):
        # runs the simulation
        print("\n\nStarting simulation")
        n_time_steps = 0
        while self.time < self.end_time:
            if (n_time_steps % self.write_freq) == 0: self.write_solution(n_time_steps)
            dt = self.calc_time_step()
            if (self.time + dt) > self.end_time:
                dt = self.end_time - self.time
            print(f"Time: {self.time:1.5e}, dt: {dt:1.5e}")
            self.do_time_step(dt)
            n_time_steps += 1
        self.write_solution(n_time_steps)
        print("Simulation done\n\n")
    
            