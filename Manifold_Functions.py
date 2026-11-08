import asset_asrl as ast
from asset_asrl.Astro.AstroModels import CR3BP
import asset_asrl.Astro.Constants as c
import numpy as np



class ManifoldSolver:
    def __init__(self, mu1=c.MuEarth, mu2=c.MuMoon, lstar=c.LD, nSeg=800):
        self.oc = ast.OptimalControl
        self.vf = ast.VectorFunctions

        self.PhaseRegion = self.oc.PhaseRegionFlags
        self.TransModes = self.oc.TranscriptionModes

        self.mu1 = mu1
        self.mu2 = mu2
        self.lstar = lstar
        self.nSeg = nSeg

        self.dt = np.pi / 10000

        self.ode = CR3BP(self.mu1, self.mu2, self.lstar)

        self.mu = self.ode.mu

        self.odeIntegrator = self.ode.integrator(self.dt)


    def solve(self, init, tf, fixedvars=[0,1,2], dS=None, dV=None, V2=None):
        tf_half = tf / 2
        
        steps = 1000
        trajguess = self.odeIntegrator.integrate_dense(init, tf_half, steps)

        odePhase = self.ode.phase(self.TransModes.LGL3)


        odePhase.setTraj(trajguess, self.nSeg)

        
        # --- STANDARD SHOOTING CONSTRAINTS ---
        for idx in fixedvars:
            odePhase.addBoundaryValue(self.PhaseRegion.Front, [idx], [trajguess[0][idx]])

        odePhase.addBoundaryValue(
                    self.PhaseRegion.Front, [3, 5, 6], [0.0, 0.0, 0.0] 
                )
        odePhase.addBoundaryValue(
                    self.PhaseRegion.Back, [1, 3, 5], [0.0, 0.0, 0.0] 
                )
            
        
        
        odePhase.optimizer.EContol = 1e-12
        odePhase.solve()


        trajSol = odePhase.returnTraj()



        return trajSol

    def monodromy(self, init, tf):


        Monodromy = self.odeIntegrator.integrate_stm(init, tf)


        return Monodromy

    def propagate(self, init, tf, steps=3000):

        traj = self.odeIntegrator.integrate_dense(init, tf, steps)


        return traj



    def generate_manifolds(self, init, tf, steps=3000, manifold_num = 20, step_size = 1e-3, multiplier = 1.0, stable=True, unstable=True, pos=True, neg=True):

        trajectories = []

        tol = 1e-5

        orbit = self.solve(init, tf)

        init_state = orbit[0]

        converged_period = orbit[-1][6] * 2

        final_state, stm = self.monodromy(init_state, converged_period)

        orbitplot = self.propagate(init_state, orbit[-1][6]*2, steps)
        trajectories.append(orbitplot)

        indices = np.linspace(0, len(orbitplot) - 1, manifold_num, endpoint=False, dtype=int)

        stm_matrix = np.array(stm)[0:6, 0:6]


        values, vectors = np.linalg.eig(stm_matrix)

        unstable_mask = (np.abs(np.imag(values)) < tol) & (np.real(values) > 1.0 + tol)
        unstable_index = np.where(unstable_mask)[0][0]

        stable_mask = (np.abs(np.imag(values)) < tol) & (np.real(values) < 1.0 - tol) & (np.real(values) > 0)
        stable_index = np.where(stable_mask)[0][0]

        v0_unstable = np.real(vectors[:, unstable_index])
        v0_stable = np.real(vectors[:, stable_index])

        integration_duration = tf * multiplier


        for idx in indices:
            state = orbitplot[idx]

            ti = state[6]

            _, stm_i = self.monodromy(init_state, ti)

            phi_i = np.array(stm_i)[0:6, 0:6]

            unstable_vector = phi_i @ v0_unstable
            stable_vector = phi_i @ v0_stable
            
            R_unstable = unstable_vector[0:3]
            Rnorm_unstable = np.linalg.norm(R_unstable)
            X_unstable = unstable_vector/Rnorm_unstable

            R_stable = stable_vector[0:3]
            Rnorm_stable = np.linalg.norm(R_stable)
            X_stable = stable_vector/Rnorm_stable

            unstable_init_pos = state.copy()
            unstable_init_pos[0:6] += step_size * X_unstable

            stable_init_pos = state.copy()
            stable_init_pos[0:6] += step_size * X_stable

            unstable_init_neg = state.copy()
            unstable_init_neg[0:6] -= step_size * X_unstable

            stable_init_neg = state.copy()
            stable_init_neg[0:6] -= step_size * X_stable

            
            
            unstable_traj_pos = self.propagate(unstable_init_pos, integration_duration, steps)
            
            unstable_traj_neg = self.propagate(unstable_init_neg, integration_duration, steps)
                
           
            
            stable_traj_pos = self.propagate(stable_init_pos, -integration_duration, steps)
            
            stable_traj_neg = self.propagate(stable_init_neg, -integration_duration, steps)
            if unstable == True:
                if pos == True:
                    trajectories.append(unstable_traj_pos)
                if neg == True:
                    trajectories.append(unstable_traj_neg)
            if stable == True:
                if pos == True:
                    trajectories.append(stable_traj_pos)
                if neg == True:
                    trajectories.append(stable_traj_neg)
                


        return trajectories