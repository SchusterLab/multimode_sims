import numpy as np
import qutip as qt
from tqdm import tqdm
import matplotlib.pyplot as plt

# There should be 2 classes , one that is just solving floquet/hamiltonian diagnalization. Let 
# call it FLoquetHamiltonian class. And the second class is its daughter class, which solves it 
# for our particular problem that involves 3 objects and couper
class FloquetHamiltonian:
    def __init__(self, H, T=None, args=None):
        """
        Initialize the FloquetHamiltonian class.

        Parameters
        ----------
        H : Qobj or list
            Time-dependent Hamiltonian or list of Hamiltonians with coefficients.
        T : float
            Floquet period (2*pi / omega).
        args : dict, optional
            Arguments for time-dependent coefficients. Must contain 'phi_AC_amp' and 'omega'.
        """
        self.H = H
        self.T = T
        if T is None:
            T = 2 * np.pi / args['omega']
        self.args = args
        self.options = qt.Options(nsteps=10000)
        self.floquet_basis = qt.FloquetBasis(H, T, args, options=self.options)

    def compute_floquet_modes(self, t=0):
        """
        Compute Floquet modes and quasienergies at a specific time.

        Parameters
        ----------
        t : float, optional
            Time at which to compute the modes. Default is 0.

        Returns
        -------
        quasienergies : np.ndarray
            Array of quasienergies.
        floquet_modes : list of Qobj
            List of Floquet mode states.
        """
        quasienergies = self.floquet_basis.e_quasi
        floquet_modes = self.floquet_basis.mode(t)
        return quasienergies, floquet_modes

    def compute_modes_over_period(self, num_points=100):
        """
        Compute Floquet modes and quasienergies over one Floquet period.

        Parameters
        ----------
        num_points : int, optional
            Number of time points to sample within one Floquet period. Default is 100.

        Returns
        -------
        energies_list : list of np.ndarray
            List of Floquet quasienergies at each time point.
        modes_list : list of Qobj
            List of Floquet modes at each time point.
        """
        times = np.linspace(0, self.T, num_points)
        energies_list = []
        modes_list = []
        for t in tqdm(times, desc="Computing Floquet modes over period"):
            energies, modes = self.compute_floquet_modes(t=t)
            energies_list.append(energies)
            modes_list.append(modes)
        return energies_list, modes_list

    def sweep_Floquet(self, args_list):
        """
        Compute Floquet modes and energies for a list of arguments.

        Parameters
        ----------
        args_list : list of dict
            List of argument dictionaries, each containing 'phi_AC_amp' and 'omega'.

        Returns
        -------
        results : dict
            Dictionary containing:
                - 'args_list': List of argument dictionaries.
                - 'energies': List of Floquet quasienergies for each set of arguments.
                - 'modes': List of Floquet modes for each set of arguments.
        """
        energies_list = []
        modes_list = []

        for args in tqdm(args_list, desc="Computing Floquet modes"):
            self.args = args  # Update args for each computation
            self.floquet_basis = qt.FloquetBasis(self.H, self.T, args, options=self.options)
            f_energies, f_modes = self.compute_floquet_modes()

            energies_list.append(f_energies)
            modes_list.append(f_modes)

        results = {
            'args_list': args_list,
            'energies': energies_list,
            'modes': modes_list
        }
        return results

    def sweep_ac_drive_amplitude(self, omega, amp_range, num_points=50):
        """
        Sweep the amplitude of the AC drive and compute Floquet modes and energies.

        Parameters
        ----------
        omega : float
            Floquet drive frequency.
        amp_range : tuple
            Range of amplitudes to sweep (min_amp, max_amp).
        num_points : int, optional
            Number of amplitude points to sample. Default is 50.

        Returns
        -------
        results : dict
            Dictionary containing:
                - 'amplitudes': List of amplitudes.
                - 'energies': List of Floquet quasienergies for each amplitude.
                - 'modes': List of Floquet modes for each amplitude.
        """
        amplitudes = np.linspace(amp_range[0], amp_range[1], num_points)
        args_list = [{'amp': amp, 'omega': omega} for amp in amplitudes]
        results = self.sweep_Floquet(args_list)
        results['amplitudes'] = amplitudes
        return results

    def sweep_ac_drive_frequency(self, amp, freq_range, num_points=50):
        """
        Sweep the frequency of the AC drive and compute Floquet modes and energies.

        Parameters
        ----------
        amp : float
            Amplitude of the AC drive.
        freq_range : tuple
            Range of frequencies to sweep (min_freq, max_freq).
        num_points : int, optional
            Number of frequency points to sample. Default is 50.

        Returns
        -------
        results : dict
            Dictionary containing:
                - 'frequencies': List of frequencies.
                - 'energies': List of Floquet quasienergies for each frequency.
                - 'modes': List of Floquet modes for each frequency.
        """
        frequencies = np.linspace(freq_range[0], freq_range[1], num_points)
        args_list = [{'amp': amp, 'omega': freq} for freq in frequencies]
        results = self.sweep_Floquet(args_list)
        results['frequencies'] = frequencies
        return results
    
    def find_max_overlap_indices(self, state, state_list, threshold=0.0, return_overlaps=False):
        """
        Given a state and a list of states, find the indices of states in the list
        with the maximum overlap above a given threshold.

        Parameters
        ----------
        state : Qobj
            The state to compare.
        state_list : list of Qobj
            List of states to compare against.
        threshold : float, optional
            Minimum overlap threshold to consider. Default is 0.0.
        return_overlaps : bool, optional
            If True, also return the filtered overlap values. Default is False.

        Returns
        -------
        indices : list of int
            List of indices of states in state_list with overlap above the threshold,
            ordered from largest to smallest overlap.
        overlaps : list of float, optional
            List of overlap values corresponding to the filtered indices, if return_overlaps is True.
        """
        overlaps = [np.abs(state.overlap(other_state)**2) for other_state in state_list]
        sorted_indices = np.argsort(overlaps)[::-1]  # Sort indices by descending overlap
        filtered_indices = [idx for idx in sorted_indices if overlaps[idx] >= threshold]
        if return_overlaps:
            filtered_overlaps = [overlaps[idx] for idx in filtered_indices]
            return filtered_indices, filtered_overlaps
        return filtered_indices
    
        




# this is analog of coupler hamiltonian but with only 2 modes and no coupler. 
# H = omega_A * a_dag * a + omega_B * b_dag * b + g_AB*sin(omega t) * (a_dag * b + a * b_dag)  
class TwoModeHamiltonian:
    def __init__(self, omega_A, omega_B, trunc):
        """
        Initialize the TwoModeHamiltonian class.

        Parameters
        ----------
        omega_A : float
            Frequency of mode A.
        omega_B : float
            Frequency of mode B.
        g_AB : float
            Coupling strength between modes A and B.
        trunc : int
            Truncation level for the Hilbert space.
        omega_drive : float
            Driving frequency.
        """
        self.omega_A = omega_A
        self.omega_B = omega_B
        # self.g_AB = g_AB
        self.trunc = trunc
        # self.omega_drive = omega_drive

        # Generate basis operators
        self.a = qt.tensor(qt.destroy(trunc), qt.qeye(trunc))
        self.a_dag = self.a.dag()
        self.b = qt.tensor(qt.qeye(trunc), qt.destroy(trunc))
        self.b_dag = self.b.dag()

        # Generate the Hamiltonian
        self.H = self.generate_H()
        
    def get_basis_state(self, n_A, n_B):
        """
        Retrieve the basis state with n_A photons in Alice and n_B photons in Bob.

        Parameters
        ----------
        n_A : int
            Number of photons in Alice's mode.
        n_B : int
            Number of photons in Bob's mode.

        Returns
        -------
        basis_state : Qobj
            The basis state with the specified photon numbers.
        """
        if n_A >= self.trunc or n_B >= self.trunc:
            raise ValueError("Photon numbers exceed truncation level.")
        basis_state = qt.tensor(
            qt.basis(self.trunc, n_A),
            qt.basis(self.trunc, n_B)
        )
        return basis_state

    def generate_H(self):
        """
        Generate the time-dependent Hamiltonian for the two-mode system.

        Returns
        -------
        H : list
            Time-dependent Hamiltonian in Floquet form.
        """
        H0 = self.omega_A * self.a_dag * self.a + self.omega_B * self.b_dag * self.b

        def coupling_coeff(t, args):
            g_AB = args['amp']
            omega_drive = args['omega']
            return g_AB * np.sin(omega_drive * t)

        H_coupling = (self.a + self.a_dag)*(self.b + self.b_dag)
        H = [H0, [H_coupling, coupling_coeff]]
        return H
    
    def get_all_basis_states(self):
        """
        Generate all computational basis states for the two-mode system.

        Returns
        -------
        basis_states : list of Qobj
            List of basis states as Qobj objects.
        ket_strings : list of str
            List of basis states as ket string representations.
        """
        basis_states = []
        ket_strings = []
        for n_A in range(self.trunc):
            for n_B in range(self.trunc):
                state = self.get_basis_state(n_A, n_B)
                basis_states.append(state)
                ket_strings.append(f"|{n_A},{n_B}>")
        return basis_states, ket_strings
    


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
        # Generate all basis states
        self.basis_states = self.generate_basis_states()
        self.generate_dressed_states()

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

    def generate_H0_with_couplings(self, lab = True):
        # Example time-dependent Hamiltonian construction for Alice and Bob
        def H_coupler_alice_coupling():
            return (self.a + self.a_dag)*(self.c + self.c_dag)

        def H_coupler_bob_coupling():
            return (self.b + self.b_dag)*(self.c + self.c_dag)

        H0 = self.generate_H0(rotating=not lab)


        if lab:
            print('In lab frame')
            # print(H0.shape)
            # print(H_coupler_alice_coupling().shape)
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

    def linc_potential_operator(self, M=3, E_J=None, E_C=None, E_L=None, rotating=False):
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
        operator = 2 * M**2 * E_J * (theta / M).cosm()

        def coeff_func(t, args):
            
            amp = args['phi_AC_amp']
            omega = args['omega']
            return np.sin(amp * np.cos(omega * t))

        return operator, coeff_func

    def static_modes(self, H):
        """
        Diagonalize the static Hamiltonian H and return the eigenvalues and eigenstates.

        Parameters
        ----------
        H : Qobj
            The Hamiltonian to diagonalize.

        Returns
        -------
        eigenvalues : np.ndarray
            Array of eigenvalues.
        eigenstates : list of Qobj
            List of eigenstate Qobjs.
        """
        eigenvalues, eigenstates = H.eigenstates()
        return np.array(eigenvalues), eigenstates
    
    # mapping between dressed states and basis states
    # Input: dressed states 
    # Output : keys are basis state indices and values are indices of dressed states
    #          and another dict where keys are dressed state indices and values are basis state indices
    
    def basis_to_dressed_mapping(self, dressed_states):
        """
        Create mappings between dressed states and basis states.

        Parameters
        ----------
        dressed_states : list of Qobj
            List of dressed states.

        Returns
        -------
        basis_to_dressed_mapping : dict
            Dictionary where keys are basis state indices and values are indices of dressed states.
        dressed_to_basis_mapping : dict
            Dictionary where keys are dressed state indices and values are basis state indices.
        """
        basis_states = self.basis_states

        basis_to_dressed_mapping = {}
        dressed_to_basis_mapping = {}
        for i, dressed_state in enumerate(dressed_states):
            overlaps = [np.abs(basis.overlap(dressed_state)**2) for basis in basis_states]
            max_idx = int(np.argmax(overlaps))
            basis_to_dressed_mapping[max_idx] = i
            dressed_to_basis_mapping[i] = max_idx

        return basis_to_dressed_mapping, dressed_to_basis_mapping
    
    # generate the dressed states 
    # Input: Nothing
    # Set self.dressed_states, self.dressed_mapping with output from above
    def generate_dressed_states(self):
        """
        Generate the dressed states by diagonalizing the static Hamiltonian
        and create a mapping between dressed states and basis states.

        Sets self.dressed_states and self.dressed_mapping.

        Returns
        -------
        dressed_states : list of Qobj
            List of dressed states.
        dressed_mapping : dict
            Mapping between basis state ket strings and dressed state indices.
        """
        self.H0 = self.generate_H0_with_couplings(lab=True)
        operator, time_func = self.linc_potential_operator()
        H_static = self.H0 + (time_func(0, args={'phi_AC_amp': 0, 'omega': 0}) * operator)

        dressed_static_energies, dressed_static_states = self.static_modes(H_static)
        basis_to_dressed_mapping_dict, dressed_to_basis_mapping_dict = self.basis_to_dressed_mapping(dressed_static_states)

        self.dressed_states = dressed_static_states
        self.basis_to_dressed_mapping_dict = basis_to_dressed_mapping_dict
        self.dressed_to_basis_mapping_dict = dressed_to_basis_mapping_dict
        self.dressed_static_energies = dressed_static_energies

        return None
    
    
    
    # Map floquet to dressed 
    # Input : floquet states, 
    # Output: floquet to dressed mapping dictionary with keys as dressed state indices and values as floquet state indices
    def dressed_to_floquet_mapping(self, floquet_states):
        """
        Create a mapping between Floquet states and dressed states.

        Parameters
        ----------
        floquet_states : list of Qobj
            List of Floquet states.

        Returns
        -------
        mapping : dict
            Dictionary where keys are dressed state indices and values are Floquet state indices.
        """
        floquet_to_dressed_mapping = {}
        dressed_to_floquet_mapping = {}
        floquet_to_dressed_mapping = {}

        for i, floquet_state in enumerate(floquet_states):
            overlaps = [np.abs(dressed_state.overlap(floquet_state)**2) for dressed_state in self.dressed_states]
            max_idx = int(np.argmax(overlaps))
            # dressed_to_floquet_mapping[max_idx] = i
            floquet_to_dressed_mapping[i] = max_idx
        # insteaf of computing maximum overlaps, iterate over every floquet state and see 
        # which dressed state it has the maximum overlap with. That is the floquet to dressed mapping
        # For dressed to floquet mapping, iterate over dressed states and see which Floquet state has the maximum overlap with it
        for i, dressed_state in enumerate(self.dressed_states):
            overlaps = [np.abs(floquet_state.overlap(dressed_state)**2) for floquet_state in floquet_states]
            max_idx = int(np.argmax(overlaps))
            dressed_to_floquet_mapping[i] = max_idx
            # floquet_to_dressed_mapping[max_idx] = i
        # for i, dressed_state in enumerate(self.dressed_states):
        #     overlaps = [np.abs(floquet_state.overlap(dressed_state)**2) for floquet_state in floquet_states]
        #     max_idx = int(np.argmax(overlaps))
        #     mapping[i] = max_idx
        return dressed_to_floquet_mapping, floquet_to_dressed_mapping
    
    # Plot floquet states in terms of dressed states 
    # input: a floquet state
    # output: a bar plot with x axis containing indices of dressed states and y axis containing the overlap with the floquet state
    def plot_floquet_state_overlap(self, floquet_state, max_total_photons=3):
        """
        Plot the overlap of a Floquet state with the dressed states, filtering based on the maximum number of photons.

        Parameters
        ----------
        floquet_state : Qobj
            The Floquet state to analyze.
        max_total_photons : int, optional
            Maximum total number of photons allowed in the mapped basis state. Default is 3.

        Returns
        -------
        None
        """
        # Compute overlaps
        overlaps = [np.abs(dressed.overlap(floquet_state)**2) for dressed in self.dressed_states]
        print("Overlaps with dressed states:", overlaps)
        plt.plot(range(len(overlaps)), overlaps, 'o', color='C0')
        # Filter dressed states based on the total photon count in the mapped basis state
        filtered_indices = []
        filtered_overlaps = []
        filtered_ket_labels = []

        for i, overlap in enumerate(overlaps):
            # Get the associated basis state index for given dressed state index
            # Count photons in the basis state
            # This is the number of photons in the dressed state, filter based on that
            basis_index = self.dressed_to_basis_mapping_dict[i]
            total_photons = self.count_total_photons(basis_index)
            if total_photons <= max_total_photons:
                filtered_indices.append(i)
                filtered_overlaps.append(overlap)
                filtered_ket_labels.append(self.basis_index_to_ket(basis_index))
       

        # Plot filtered overlaps
        plt.figure(figsize=(10, 6))
        plt.bar(filtered_indices, filtered_overlaps, color='C0')
        plt.xlabel('Dressed State Index')
        plt.ylabel('Overlap Probability')
        plt.title('Overlap of Floquet State with Dressed States (Filtered by Total Photons)')
        plt.xticks(filtered_indices, [f"~{label}" for label in filtered_ket_labels], rotation=90)
        plt.tight_layout()
        plt.show()
    
    # Plot a dressed state in terms of the basis states
    # Input: a dressed state 
    # Output: a bar plot with x axis containing indices of basis states and y axis containing the overlap with the dressed state
    def count_total_photons(self, basis_state_index):
        """
        Given a basis state index, return the total number of photons across all modes.

        Parameters
        ----------
        basis_state_index : int
            Index of the basis state.

        Returns
        -------
        total_photons : int
            Total number of photons in the basis state.
        """
        n_A = basis_state_index // (self.trunc ** 2)
        n_B = (basis_state_index // self.trunc) % self.trunc
        n_C = basis_state_index % self.trunc
        return n_A + n_B + n_C

    def plot_dressed_state_overlap(self, dressed_state, max_total_photons=3):
        """
        Plot the overlap of a dressed state with the basis states, ignoring states with more than max_total_photons.

        Parameters
        ----------
        dressed_state : Qobj
            The dressed state to analyze.
        max_total_photons : int, optional
            Maximum total number of photons allowed in a basis state. Default is 3.

        Returns
        -------
        None
        """
        # Compute overlaps
        overlaps = [np.abs(basis.overlap(dressed_state)**2) for basis in self.basis_states]

        # Filter basis states based on total photon count
        filtered_indices = []
        filtered_overlaps = []
        for i, overlap in enumerate(overlaps):
            if self.count_total_photons(i) <= max_total_photons:
                filtered_indices.append(i)
                filtered_overlaps.append(overlap)

        # Plot filtered overlaps
        plt.figure(figsize=(10, 6))
        plt.bar(filtered_indices, filtered_overlaps, color='C0')
        plt.xlabel('Basis State Index')
        plt.ylabel('Overlap Probability')
        plt.title('Overlap of Dressed State with Basis States (Filtered by Total Photons)')
        plt.xticks(filtered_indices, [self.basis_index_to_ket(i) for i in filtered_indices], rotation=90)
        plt.tight_layout()
        plt.show()

    def generate_basis_states(self):
        """
        Generate all computational basis states for the system.

        Returns
        -------
        basis_states : list of Qobj
            List of basis states as Qobj objects.
        """
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
        return basis_states

    
    def floquet_modes(self,  H, omega, args=None, rotating=False, t = None):
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
        
        T = 2 * np.pi / omega

        options = qt.Options(nsteps=10000)
        floquet_basis = qt.FloquetBasis(H, T, args, options=options)
        f_energies = floquet_basis.e_quasi
        if t is None:
            t = T/2
        f_modes = floquet_basis.mode(t)

        # f_modes_0, f_energies = qt.floquet_modes(H, T, args)
        # # at later time t 
        # # if t is None:
        # #     t = 0
        # t = T/2
        # modes = qt.floquet_modes_t(f_modes_0, f_energies, t, H, T, args)
        return f_energies, f_modes
    
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

    # 
    # def compute_floquet_and_mappings(self, H, args_list):
    #     """
    #     Compute Floquet modes, energies, and mappings for a list of arguments.

    #     Parameters
    #     ----------
    #     H : list
    #         Time-dependent Hamiltonian in Floquet form.
    #     args_list : list of dict
    #         List of argument dictionaries, each containing 'phi_AC_amp' and 'omega'.

    #     Returns
    #     -------
    #     results : dict
    #         Dictionary containing:
    #             - 'args_list': List of argument dictionaries.
    #             - 'energies': List of Floquet quasienergies for each set of arguments.
    #             - 'modes': List of Floquet modes for each set of arguments.
    #             - 'mappings': List of mappings (dressed to Floquet and Floquet to dressed) for each set of arguments.
    #     """
    #     energies_list = []
    #     modes_list = []
    #     dressed_to_floquet_mapping_list = []
    #     floquet_to_dressed_mapping_list = []

    #     for args in tqdm(args_list, desc="Computing Floquet modes and mappings"):
    #         f_energies, f_modes = self.floquet_modes(H, args['omega'], args=args)
    #         dressed_to_floquet_mapping, floquet_to_dressed_mapping = self.dressed_to_floquet_mapping(f_modes)

    #         energies_list.append(f_energies)
    #         modes_list.append(f_modes)
    #         dressed_to_floquet_mapping_list.append(dressed_to_floquet_mapping)
    #         floquet_to_dressed_mapping_list.append(floquet_to_dressed_mapping)

    #     results = {
    #         'args_list': args_list,
    #         'energies': energies_list,
    #         'modes': modes_list,
    #         'mappings': {
    #             'dressed_to_floquet': dressed_to_floquet_mapping_list,
    #             'floquet_to_dressed': floquet_to_dressed_mapping_list
    #         }
    #     }
    #     return results

    # def sweep_ac_drive_amplitude(self, H, omega, amp_range, num_points=50):
    #     """
    #     Sweep the amplitude of the AC drive and compute Floquet modes, energies, and mappings.

    #     Parameters
    #     ----------
    #     H : list
    #         Time-dependent Hamiltonian in Floquet form.
    #     omega : float
    #         Floquet drive frequency.
    #     amp_range : tuple
    #         Range of amplitudes to sweep (min_amp, max_amp).
    #     num_points : int, optional
    #         Number of amplitude points to sample. Default is 50.

    #     Returns
    #     -------
    #     results : dict
    #         Dictionary containing:
    #             - 'amplitudes': List of amplitudes.
    #             - 'energies': List of Floquet quasienergies for each amplitude.
    #             - 'modes': List of Floquet modes for each amplitude.
    #             - 'mappings': List of mappings (dressed to Floquet and Floquet to dressed) for each amplitude.
    #     """
    #     amplitudes = np.linspace(amp_range[0], amp_range[1], num_points)
    #     args_list = [{'phi_AC_amp': amp, 'omega': omega} for amp in amplitudes]
    #     results = self.compute_floquet_and_mappings(H, args_list)
    #     results['amplitudes'] = amplitudes
    #     return results

    # def sweep_ac_drive_frequency(self, H, amp, freq_range, num_points=50):
    #     """
    #     Sweep the frequency of the AC drive and compute Floquet modes, energies, and mappings.

    #     Parameters
    #     ----------
    #     H : list
    #         Time-dependent Hamiltonian in Floquet form.
    #     amp : float
    #         Amplitude of the AC drive.
    #     freq_range : tuple
    #         Range of frequencies to sweep (min_freq, max_freq).
    #     num_points : int, optional
    #         Number of frequency points to sample. Default is 50.

    #     Returns
    #     -------
    #     results : dict
    #         Dictionary containing:
    #             - 'frequencies': List of frequencies.
    #             - 'energies': List of Floquet quasienergies for each frequency.
    #             - 'modes': List of Floquet modes for each frequency.
    #             - 'mappings': List of mappings (dressed to Floquet and Floquet to dressed) for each frequency.
    #     """
    #     frequencies = np.linspace(freq_range[0], freq_range[1], num_points)
    #     args_list = [{'phi_AC_amp': amp, 'omega': freq} for freq in frequencies]
    #     results = self.compute_floquet_and_mappings(H, args_list)
    #     results['frequencies'] = frequencies
    #     return results
    # # Sweep the amplitude of the AC drive and compute the Floquet modes and energies.
    # # for each mode, also store the mapping from dressed states to Floquet states and vice versa
    # # return the amplitudes, energies, modes, and the mappings in a dictionary
    # def sweep_ac_drive_amplitude(self, H, omega, amp_range, num_points=50):
    #     """
    #     Sweep the amplitude of the AC drive and compute Floquet modes, energies, and mappings.

    #     Parameters
    #     ----------
    #     H : list
    #         Time-dependent Hamiltonian in Floquet form.
    #     omega : float
    #         Floquet drive frequency.
    #     amp_range : tuple
    #         Range of amplitudes to sweep (min_amp, max_amp).
    #     num_points : int, optional
    #         Number of amplitude points to sample. Default is 50.

    #     Returns
    #     -------
    #     results : dict
    #         Dictionary containing:
    #             - 'amplitudes': List of amplitudes.
    #             - 'energies': List of Floquet quasienergies for each amplitude.
    #             - 'modes': List of Floquet modes for each amplitude.
    #             - 'mappings': List of mappings (dressed to Floquet and Floquet to dressed) for each amplitude.
    #     """
    #     amplitudes = np.linspace(amp_range[0], amp_range[1], num_points)
    #     energies_list = []
    #     modes_list = []
    #     dressed_to_floquet_mapping_list = []
    #     floquet_to_dressed_mapping_list = []

    #     for amp in tqdm(amplitudes, desc="Sweeping AC drive amplitude"):
    #         args = {'phi_AC_amp': amp, 'omega': omega}
    #         f_energies, f_modes = self.floquet_modes(H, omega, args=args)
    #         dressed_to_floquet_mapping, floquet_to_dressed_mapping = self.dressed_to_floquet_mapping(f_modes)

    #         energies_list.append(f_energies)
    #         modes_list.append(f_modes)
    #         dressed_to_floquet_mapping_list.append(dressed_to_floquet_mapping)
    #         floquet_to_dressed_mapping_list.append(floquet_to_dressed_mapping)

    #     results = {
    #         'amplitudes': amplitudes,
    #         'energies': energies_list,
    #         'modes': modes_list,
    #         'mappings': {
    #             'dressed_to_floquet': dressed_to_floquet_mapping_list,
    #             'floquet_to_dressed': floquet_to_dressed_mapping_list
    #         }
    #     }
    #     return results

    # def plot_floquet_energies_over_period(self, energies_list, modes_list, max_photons_per_mode=2, ax=None):
    #     """
    #     Plot the Floquet quasienergies for each basis state over time.

    #     Parameters
    #     ----------
    #     en_for_state : np.ndarray
    #         2D array of shape (num_states, num_points) with energies for each state over time.
    #     max_photons_per_mode : int, optional
    #         Maximum number of photons allowed per mode for a state to be plotted.
    #     ax : matplotlib.axes.Axes, optional
    #         Axis to plot on. If None, creates a new figure and axis.
    #     """
    #     import matplotlib.pyplot as plt

    #     num_points, num_states = np.shape(energies_list)
        
    #     if ax is None:
    #         fig, ax = plt.subplots()
        
    #     big_mapping_list = [self.floquet_to_basis_mapping(modes) for modes in modes_list]
    #     state_idxs = list(big_mapping_list[0].keys())
    #     en_for_state = np.zeros(( num_states, num_points))
    #     for j, key in enumerate(state_idxs):
    #         print(big_mapping_list[0][key])
    #         en_for_state[j, :] = [energies_list[i][big_mapping_list[i][key]] for i in range(num_points)]
        
    #     # kets = self.convert_basis_index_list_to_ket_list(big_mapping_list[0].keys())

    #     for j in state_idxs:
    #         n_A = j // (self.trunc ** 2)
    #         n_B = (j // self.trunc) % self.trunc
    #         n_C = j % self.trunc
    #         if max_photons_per_mode is not None and (n_A > max_photons_per_mode or n_B > max_photons_per_mode or n_C > max_photons_per_mode):
    #             continue
    #         label = self.basis_index_to_ket(j)
    #         ax.plot(range(num_points), en_for_state[j, :], label=label)

    #     ax.set_xlabel('Time Index')
    #     ax.set_ylabel('Quasienergy')
    #     ax.set_title('Floquet Quasienergies Over Time')
    #     ax.legend(fontsize=8, loc='best', ncol=2)
    #     plt.tight_layout()
    #     plt.show()
    
    # def plot_floquet_energies(self, floquet_states, quasienergies, ax=None, max_photons_per_mode=2):
    #     """
    #     Plot quasienergies of Floquet states with text labels indicating their basis states.

    #     Parameters
    #     ----------
    #     floquet_states : list of Qobj
    #         List of Floquet mode states (Qobj) to plot.
    #     quasienergies : list or np.ndarray
    #         Quasienergies corresponding to each Floquet state.
    #     ax : matplotlib.axes.Axes, optional
    #         Axis to plot on. If None, creates a new figure and axis.
    #     max_photons_per_mode : int or None, optional
    #         If set, ignore states where any mode has more than this number of photons.
    #     """
    #     import matplotlib.pyplot as plt

    #     if ax is None:
    #         fig, ax = plt.subplots()

    #     # Use the mapping from Floquet state index to basis state index
    #     mapping = self.floquet_to_basis_mapping(floquet_states)

    #     for i, energy in enumerate(quasienergies):
    #         basis_idx = mapping[i]
    #         n_A = basis_idx // (self.trunc ** 2)
    #         n_B = (basis_idx // self.trunc) % self.trunc
    #         n_C = basis_idx % self.trunc
    #         if max_photons_per_mode is not None and (n_A > max_photons_per_mode or n_B > max_photons_per_mode or n_C > max_photons_per_mode):
    #             continue
    #         label = self.basis_index_to_ket(basis_idx)
    #         ax.plot(i, energy, 'o', color='C0')
    #         ax.text(i, energy, label, ha='center', va='bottom', fontsize=9, rotation=90)

    #     ax.set_xlabel('Floquet State Index')
    #     ax.set_ylabel('Quasienergy')
    #     ax.set_title('Floquet Quasienergies with Basis Labels')
    #     plt.tight_layout()
    #     plt.show()


    # def floquet_to_basis_mapping(self, floquet_states):
    #     """
    #     Create a dictionary mapping Floquet states to original basis states
    #     based on the largest overlap.

    #     Parameters
    #     ----------
    #     floquet_states : list of Qobj
    #         List of Floquet mode states.

    #     Returns
    #     -------
    #     mapping : dict
    #         Dictionary mapping Floquet state index to basis state index.
    #     """
    #     # Construct computational basis states
    #     dim = self.trunc
    #     basis_states = []
    #     for dim_a in range(self.trunc):
    #         for dim_b in range(self.trunc):
    #             for dim_c in range(self.trunc):
    #                 state = qt.tensor(
    #                     qt.basis(dim, dim_a),
    #                     qt.basis(dim, dim_b),
    #                     qt.basis(dim, dim_c)
    #                 )
    #                 basis_states.append(state)
    #     mapping = {}
    #     for i, floquet_state in enumerate(floquet_states):
    #         overlaps = [np.abs(basis.overlap(floquet_state)**2) for basis in basis_states]
    #         max_idx = int(np.argmax(overlaps))
    #         mapping[i] = max_idx
    #     return mapping

    # given a state and list of states, find the index of state with maximum overlap, 
    # if there are multiple above a threshold say 0.9, returna list of indices ordered from largest to smallest overlap
    def find_max_overlap_indices(self, state, state_list, threshold=0.0):
        """
        Given a state and a list of states, find the indices of states in the list
        with the maximum overlap above a given threshold.

        Parameters
        ----------
        state : Qobj
            The state to compare.
        state_list : list of Qobj
            List of states to compare against.
        threshold : float, optional
            Minimum overlap threshold to consider. Default is 0.9.

        Returns
        -------
        indices : list of int
            List of indices of states in state_list with overlap above the threshold,
            ordered from largest to smallest overlap.
        """
        overlaps = [np.abs(state.overlap(other_state)**2) for other_state in state_list]
        sorted_indices = np.argsort(overlaps)[::-1]  # Sort indices by descending overlap
        filtered_indices = [idx for idx in sorted_indices if overlaps[idx] >= threshold]
        return filtered_indices
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
   
