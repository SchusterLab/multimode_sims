import numpy as np
import qutip as qt
from tqdm import tqdm



class CouplerHamiltonian:
    
    def __init__(self, trunc, omega_A, omega_B, omega_C, g_AC, g_BC, E_C, E_L, E_J):
        print("Initializing CouplerHamiltonian with truncation:", trunc)
        self.trunc = trunc
        self.omega_A = omega_A
        self.omega_B = omega_B
        self.omega_C = omega_C
        self.g_AC = g_AC
        self.g_BC = g_BC
        self.E_C = E_C
        self.E_L = E_L
        self.E_J = E_J
        self.a, self.a_dag, self.b, self.b_dag, self.c, self.c_dag = self.generate_basis_operators()

    def generate_basis_operators(self):
        I = qt.qeye(self.trunc)
        a = qt.tensor(qt.destroy(self.trunc), I, I)
        a_dag = a.dag()
        b = qt.tensor(I, qt.destroy(self.trunc), I)
        b_dag = b.dag()
        c = qt.tensor(I, I, qt.destroy(self.trunc))
        c_dag = c.dag()
        return a, a_dag, b, b_dag, c, c_dag

    def generate_H0(self, rotating=False):
        if rotating:
            return qt.qeye(self.trunc ** 3)
            
        
        H = (
            self.omega_A * self.a_dag * self.a
            + self.omega_B * self.b_dag * self.b
            + self.omega_C * self.c_dag * self.c
            + self.g_AC * (self.a_dag * self.c + self.a * self.c_dag)
            + self.g_BC * (self.b_dag * self.c + self.b * self.c_dag)
        )
        return H

    def generate_H0_with_couplings(self, rotating=False):
        # Example time-dependent Hamiltonian construction for Alice and Bob
        def H_coupler_alice_coupling():
            return (self.a + self.a_dag)*(self.c + self.c_dag)

        def H_coupler_bob_coupling():
            return (self.b + self.b_dag)*(self.c + self.c_dag)

        H0 = self.generate_H0(rotating=rotating)

        if rotating:
            # In rotating frame, just add static g_AC coupling
            print(H0.shape)
            print(H_coupler_alice_coupling().shape)
            H = H0 +self.g_AC*H_coupler_alice_coupling() + self.g_BC*H_coupler_bob_coupling()
            return H
        else:
            alice_coupling = [
                [H_coupler_alice_coupling(), lambda t, args: self.g_AC * np.exp(1j * (self.omega_A - self.omega_C) * t)],
                [H_coupler_alice_coupling().dag(), lambda t, args: self.g_AC * np.exp(-1j * (self.omega_A - self.omega_C) * t)]
            ]
            bob_coupling = [
                [H_coupler_bob_coupling(), lambda t, args: self.g_BC * np.exp(1j * (self.omega_B - self.omega_C) * t)],
                [H_coupler_bob_coupling().dag(), lambda t, args: self.g_BC * np.exp(-1j * (self.omega_B - self.omega_C) * t)]
            ]
            H = [H0] + alice_coupling + bob_coupling
            return H


    def get_theta_operator(self, E_C=None, E_L=None):
        """
        Construct the theta operator for the coupler mode:
            theta = theta_zpt * (c + c†)
            where theta_zpt = (2 * E_C / E_L) ** 0.25

        Parameters
        ----------
        E_C : float, optional
            Charging energy. If None, uses self.E_C.
        E_L : float, optional
            Inductive energy. If None, uses self.E_L.

        Returns
        -------
        theta : Qobj
            The theta operator.
        """
        if E_C is None:
            E_C = self.E_C
        if E_L is None:
            E_L = self.E_L
        theta_zpt = (2 * E_C / E_L) ** 0.25
        return theta_zpt * (self.c + self.c_dag)
    
    # Modified linc_potential_operator to return the operator aspect and a cos(omega*t) function

    def linc_potential_operator(self, omega, phi_AC_amp=0.1, M=3, E_J=None, E_C=None, E_L=None, rotating=False):
        """
        Generate the LINC potential energy operator as a time-dependent Hamiltonian:
            U = 2 * M**2 * E_J * phi_AC_amp * cos(omega * t) * cos(theta / M)
        where theta = theta_zpt * (c + c†), theta_zpt = (2 * E_C / E_L) ** 0.25

        Parameters
        ----------
        M : int
            Integer parameter related to the number of junctions or modes.
        E_J : float
            Josephson energy.
        phi_AC_amp : float
            Amplitude of the phase drive.
        omega : float
            Frequency of the drive.
        E_C : float, optional
            Charging energy. If None, uses self.E_C.
        E_L : float, optional
            Inductive energy. If None, uses self.E_L.
        rotating : bool, optional
            If True, indicates rotating frame (not implemented).

        Returns
        -------
        operator : Qobj
            The operator part (cos(theta / M)).
        coeff_func : function
            Function of t and args returning 2 * M**2 * E_J * phi_AC_amp * cos(omega * t).
        """
        if rotating:
            raise NotImplementedError("linc_potential_operator in rotating frame is not implemented yet.")
        if E_C is None:
            E_C = self.E_C
        if E_L is None:
            E_L = self.E_L
        if E_J is None:
            E_J = self.E_J
        theta = self.get_theta_operator(E_C, E_L)
        operator = 2 * M**2 * E_J * phi_AC_amp * (theta / M).cosm()

        def coeff_func(t, args):
            return np.cos(omega * t)

        return operator, coeff_func
    
    def floquet_modes(self, omega, H, args=None, rotating=False, t = None):
        """
        Compute Floquet modes and quasienergies for the time-dependent Hamiltonian.
        The period T is set to the minimum of the inverse frequency differences.

        Parameters
        ----------
        args : dict, optional
            Arguments for time-dependent coefficients.
        rotating : bool, optional
            Whether to use the rotating frame Hamiltonian.

        Returns
        -------
        quasienergies : np.ndarray
            Array of quasienergies.
        floquet_modes : list of Qobj
            List of Floquet mode states.
        """
        # Compute frequency differences
        freq_diffs = [
            abs(self.omega_A - self.omega_C),
            abs(self.omega_B - self.omega_C),
            abs(self.omega_A - self.omega_B)
        ]
        # Avoid division by zero
        # freq_diffs = [fd for fd in freq_diffs if fd > 0]
        # if not freq_diffs:
        #     raise ValueError("All frequency differences are zero; cannot determine Floquet period.")
        # min_freq_diff = min(freq_diffs)
        T = 2 * np.pi / omega

        options = qt.Options(nsteps=10000)
        floquet_basis = qt.FloquetBasis(H, T, args, options=options)
        f_energies = floquet_basis.e_quasi
        if t is None:
            t = T/2
        f_modes_0 = floquet_basis.mode(t)
        return f_energies, f_modes_0
    
    def compute_states_over_floquet_period(self, omega, H, num_points=100):
        """
        Computes the Floquet states and their corresponding energies over one Floquet period.
        This method evaluates the Floquet modes and energies at evenly spaced time points
        within a single period of the driving frequency. It also computes the mapping from
        Floquet modes to the original basis at each time point.
        Args:
            omega (float): The driving frequency.
            H (callable or np.ndarray): The system Hamiltonian or a function returning the Hamiltonian at time t.
            num_points (int, optional): Number of time points to sample within one Floquet period. Default is 100.
        Returns:
            tuple:
                energies_list (list of np.ndarray): List of Floquet quasi-energies at each time point.
                modes_list (list of np.ndarray): List of Floquet modes at each time point.
                big_mapping_list (list): List of mappings from Floquet modes to the original basis at each time point.
        """
        

        T = 2 * np.pi / omega
        energies_list = []
        modes_list = []
        for t in tqdm(np.linspace(0, T, num_points), desc="Computing Floquet states"):
            energies, modes = self.floquet_modes(omega, H, t=t)
            energies_list.append(energies)
            modes_list.append(modes)
        
        # get mappings for all the modes
        # big_mapping_list = [self.floquet_to_basis_mapping(modes) for modes in modes_list]
        return energies_list, modes_list

    def plot_floquet_energies_over_period(self, energies_list, modes_list, max_photons_per_mode=2, ax=None):
        """
        Plot the Floquet quasienergies for each basis state over time.

        Parameters
        ----------
        en_for_state : np.ndarray
            2D array of shape (num_states, num_points) with energies for each state over time.
        max_photons_per_mode : int, optional
            Maximum number of photons allowed per mode for a state to be plotted.
        ax : matplotlib.axes.Axes, optional
            Axis to plot on. If None, creates a new figure and axis.
        """
        import matplotlib.pyplot as plt

        num_points, num_states = np.shape(energies_list)
        
        if ax is None:
            fig, ax = plt.subplots()
        
        big_mapping_list = [self.floquet_to_basis_mapping(modes) for modes in modes_list]
        state_idxs = list(big_mapping_list[0].keys())
        en_for_state = np.zeros(( num_states, num_points))
        for j, key in enumerate(state_idxs):
            print(big_mapping_list[0][key])
            en_for_state[j, :] = [energies_list[i][big_mapping_list[i][key]] for i in range(num_points)]
        
        # kets = self.convert_basis_index_list_to_ket_list(big_mapping_list[0].keys())

        for j in state_idxs:
            n_A = j // (self.trunc ** 2)
            n_B = (j // self.trunc) % self.trunc
            n_C = j % self.trunc
            if max_photons_per_mode is not None and (n_A > max_photons_per_mode or n_B > max_photons_per_mode or n_C > max_photons_per_mode):
                continue
            label = self.basis_index_to_ket(j)
            ax.plot(range(num_points), en_for_state[j, :], label=label)

        ax.set_xlabel('Time Index')
        ax.set_ylabel('Quasienergy')
        ax.set_title('Floquet Quasienergies Over Time')
        ax.legend(fontsize=8, loc='best', ncol=2)
        plt.tight_layout()
        plt.show()
    
    def plot_floquet_energies(self, floquet_states, quasienergies, ax=None, max_photons_per_mode=2):
        """
        Plot quasienergies of Floquet states with text labels indicating their basis states.

        Parameters
        ----------
        floquet_states : list of Qobj
            List of Floquet mode states (Qobj) to plot.
        quasienergies : list or np.ndarray
            Quasienergies corresponding to each Floquet state.
        ax : matplotlib.axes.Axes, optional
            Axis to plot on. If None, creates a new figure and axis.
        max_photons_per_mode : int or None, optional
            If set, ignore states where any mode has more than this number of photons.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()

        # Use the mapping from Floquet state index to basis state index
        mapping = self.floquet_to_basis_mapping(floquet_states)

        for i, energy in enumerate(quasienergies):
            basis_idx = mapping[i]
            n_A = basis_idx // (self.trunc ** 2)
            n_B = (basis_idx // self.trunc) % self.trunc
            n_C = basis_idx % self.trunc
            if max_photons_per_mode is not None and (n_A > max_photons_per_mode or n_B > max_photons_per_mode or n_C > max_photons_per_mode):
                continue
            label = self.basis_index_to_ket(basis_idx)
            ax.plot(i, energy, 'o', color='C0')
            ax.text(i, energy, label, ha='center', va='bottom', fontsize=9, rotation=90)

        ax.set_xlabel('Floquet State Index')
        ax.set_ylabel('Quasienergy')
        ax.set_title('Floquet Quasienergies with Basis Labels')
        plt.tight_layout()
        plt.show()


    def floquet_to_basis_mapping(self, floquet_states):
        """
        Create a dictionary mapping Floquet states to original basis states
        based on the largest overlap.

        Parameters
        ----------
        floquet_states : list of Qobj
            List of Floquet mode states.

        Returns
        -------
        mapping : dict
            Dictionary mapping Floquet state index to basis state index.
        """
        # Construct computational basis states
        dim = self.trunc
        basis_states = []
        for dim_a in range(self.trunc):
            for dim_b in range(self.trunc):
                for dim_c in range(self.trunc):
                    state = qt.tensor(
                        qt.basis(dim, dim_a),
                        qt.basis(dim, dim_b),
                        qt.basis(dim, dim_c)
                    )
                    basis_states.append(state)
        mapping = {}
        for i, floquet_state in enumerate(floquet_states):
            overlaps = [np.abs(basis.overlap(floquet_state)**2) for basis in basis_states]
            max_idx = int(np.argmax(overlaps))
            mapping[i] = max_idx
        return mapping

    def convert_basis_index_list_to_ket_list(self, indices):
        """
        Given a list of basis state indices, return a list of their ket string representations.

        Parameters
        ----------
        indices : list of int
            List of basis state indices.

        Returns
        -------
        ket_list : list of str
            List of ket string representations.
        """
        return [self.basis_index_to_ket(idx) for idx in indices]
    def basis_index_to_ket(self, index):
        """
        Given a basis state index, return the ket string representation |n_A, n_B, n_C>.

        Parameters
        ----------
        index : int
            Basis state index.

        Returns
        -------
        ket_str : str
            String representation, e.g., '|0,1,2>'.
        """
        n_A = index // (self.trunc ** 2)
        n_B = (index // self.trunc) % self.trunc
        n_C = index % self.trunc
        return f"|{n_A},{n_B},{n_C}>"

    def ket_to_basis_index(self, ket_str):
        """
        Given a ket string representation, return the basis state index.

        Parameters
        ----------
        ket_str : str
            String representation, e.g., '|0,1,2>'.

        Returns
        -------
        index : int
            Basis state index.
        """
        # Remove '|' and '>' and split by ','
        nums = ket_str.strip('|>').split(',')
        n_A, n_B, n_C = map(int, nums)
        return n_A * (self.trunc ** 2) + n_B * self.trunc + n_C
   
