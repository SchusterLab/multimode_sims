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
        # self.floquet_basis = qt.FloquetBasis(H, T, args, options=self.options)

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
        # self.floquet_basis = qt.FloquetBasis(H, T, args, options=self.options)
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
            T = 2 * np.pi / args['omega'] if 'omega' in args else self.T
            self.floquet_basis = qt.FloquetBasis(self.H, T, args, options=self.options)
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
    
    # given a list of modes, one list for each swept parameter point, find the max overlap of target mode with the list . 
    # then plot the overlaps as a function of the swept parameter
    def plot_overlaps_states(self, swept_list, target_states, target_state_labels, modes_swept_list, 
                             xlabel="Swept Parameter", ylabel='Overlap', find_max=True, find_min=False):
        """
        Plot overlaps for multiple target states, each with separate plots.

        Parameters
        ----------
        swept_list : list
            List of swept parameter values.
        target_states : list of Qobj
            List of target states to compute overlaps for.
        target_state_labels : list of str
            List of labels for the target states.
        modes_swept_list : list of list of Qobj
            List of lists of modes corresponding to each swept parameter value.
        xlabel : str, optional
            Label for the x-axis. Default is "Swept Parameter".
        ylabel : str, optional
            Label for the y-axis. Default is "Overlap".
        find_max : bool, optional
            If True, find and mark the maximum overlap. Default is True.
        find_min : bool, optional
            If True, find and mark the minimum overlap. Default is False.

        Returns
        -------
        results : list of dict
            List of dictionaries containing swept parameters, overlaps, indices, and max/min indices for each target state.
        """
        results = []

        for target_state, target_label in zip(target_states, target_state_labels):
            overlap_idxs = []
            overlap_values = []
            for i in range(len(swept_list)):
                alice_overlap_indices, alice_overlap = self.find_max_overlap_indices(target_state, modes_swept_list[i], return_overlaps=True)
                overlap_idxs.append(alice_overlap_indices[0])  # Take the first index (max overlap)
                overlap_values.append(alice_overlap[0])

            # plt.figure(figsize=(10, 6))
            fig, ax = plt.subplots(2, 1, figsize=(10, 12))

            # Plot overlaps
            ax[0].plot(swept_list, overlap_values, marker='o', linestyle='-', color='C0')
            ax[0].set_xlabel(xlabel)
            ax[0].set_ylabel(ylabel)
            ax[0].set_title(f'Max Overlap of Target State ({target_label}) with Swept Modes')
            ax[0].grid()

            # Plot state indices
            ax[1].plot(swept_list, overlap_idxs, marker='o', linestyle='-', color='C1')
            ax[1].set_xlabel(xlabel)
            ax[1].set_ylabel('State Index')
            ax[1].set_title(f'State Index of Max Overlap for Target State ({target_label})')
            ax[1].grid()

            # Highlight max and min overlap points
            max_overlap_idx = None
            min_overlap_idx = None

            if find_max:
                max_overlap_idx = np.argmax(overlap_values)
                print(f"Max overlap at index {max_overlap_idx}: {overlap_values[max_overlap_idx]} for swept parameter {swept_list[max_overlap_idx]}")
                ax[0].axvline(x=swept_list[max_overlap_idx], color='r', linestyle='--', label='Max Overlap')
            if find_min:
                min_overlap_idx = np.argmin(overlap_values)
                print(f"Min overlap at index {min_overlap_idx}: {overlap_values[min_overlap_idx]} for swept parameter {swept_list[min_overlap_idx]}")
                ax[0].axvline(x=swept_list[min_overlap_idx], color='g', linestyle='--', label='Min Overlap')

            ax[0].legend()
            # plt.tight_layout()
            # plt.show()
            # plt.plot(swept_list, overlap_values, marker='o', linestyle='-', color='C0')

            # plt.xlabel(xlabel)
            # plt.ylabel(ylabel)
            # plt.title(f'Max Overlap of Target State ({target_label}) with Swept Modes')
            # plt.grid()

            # max_overlap_idx = None
            # min_overlap_idx = None

            # if find_max:
            #     max_overlap_idx = np.argmax(overlap_values)
            #     print(f"Max overlap at index {max_overlap_idx}: {overlap_values[max_overlap_idx]} for swept parameter {swept_list[max_overlap_idx]}")
            #     plt.axvline(x=swept_list[max_overlap_idx], color='r', linestyle='--', label='Max Overlap')
            # if find_min:
            #     min_overlap_idx = np.argmin(overlap_values)
            #     print(f"Min overlap at index {min_overlap_idx}: {overlap_values[min_overlap_idx]} for swept parameter {swept_list[min_overlap_idx]}")
            #     plt.axvline(x=swept_list[min_overlap_idx], color='g', linestyle='--', label='Min Overlap')

            # plt.legend()
            # plt.tight_layout()
            # plt.show()

            results.append({
                'swept_list': swept_list,
                'overlap_values': overlap_values,
                'overlap_idxs': overlap_idxs,
                'max_overlap_idx': max_overlap_idx,
                'min_overlap_idx': min_overlap_idx
            })

        return results
           




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

    def generate_H0(self, rotating=False, no_coupling=False):
        if rotating:
            return qt.qeye(self.trunc ** 3)
        
        H = (
            self.omega_A * self.a_dag * self.a
            + self.omega_B * self.b_dag * self.b
            + self.omega_C * self.c_dag * self.c
        )
        if no_coupling:
            print("No couplings added to Hamiltonian")
            
        if not no_coupling:
            print("Adding couplings to Hamiltonian")
            H += (
                self.g_AC * (self.a_dag - self.a) * (self.c_dag - self.c)
                + self.g_BC * (self.b_dag - self.b) * (self.c_dag - self.c)
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
            
            amp = args['amp']
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
        H_static = self.H0 + (time_func(0, args={'amp': 0, 'omega': 0}) * operator)

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
    # ---------------------Finding effective g between states---------------------
    # Compute Floquet modes and energies, use either dressed states to idenitfy alice bob or floquet driven uncoupled states 
    def calculate_coupled_uncoupled_dressed_overlaps_and_energies(self, floquet_class_uncoup, floquet_class_coup,  return_args_uncoup, return_args_coup, state_label):
        """
        Calculates the overlaps and energies between a specified state in the coupled system and corresponding states in both uncoupled and dressed bases across a set of frequencies.
        For each frequency, this method:
            - Finds the overlap and energy of the coupled state with the closest uncoupled state.
            - Finds the overlap and energy of the coupled state with the closest dressed state.
        Args:
            floquet_class_uncoup: An instance representing the uncoupled Floquet system, providing methods to find overlaps.
            floquet_class_coup: An instance representing the coupled Floquet system, providing methods to find overlaps.
            return_args_uncoup (dict): Dictionary containing results for the uncoupled system, including 'modes' and other relevant data.
            return_args_coup (dict): Dictionary containing results for the coupled system, including 'modes', 'energies', and 'frequencies'.
            state_label: Label or identifier for the state of interest in the coupled system.
        Returns:
            tuple: A tuple containing four lists:
                - overlaps_coup_with_uncoupled (list): Overlaps of the coupled state with the closest uncoupled states for each frequency.
                - energies_coup_with_uncoupled (list): Energies of the closest uncoupled states for each frequency.
                - overlaps_coup_with_dressed (list): Overlaps of the coupled state with the closest dressed states for each frequency.
                - energies_coup_with_dressed (list): Energies of the closest dressed states for each frequency.
        """
        
        overlaps_coup_with_uncoupled = []
        energies_coup_with_uncoupled = []
        overlaps_coup_with_dressed = []
        energies_coup_with_dressed = []
        coup_ham = self

        for idx, freq in enumerate(return_args_coup['frequencies']):
            state = coup_ham.basis_states[coup_ham.ket_to_basis_index(state_label)]

            # Find the state in the uncoupled modes
            overlap_idxs, overlaps = floquet_class_uncoup.find_max_overlap_indices(
                state, return_args_uncoup['modes'][idx], return_overlaps=True
            )
            driven_uncoupled_state = return_args_uncoup['modes'][idx][overlap_idxs[0]]

            # Find the states in driven coupled modes closest to the uncoupled states
            overlap_idxs, overlaps = floquet_class_coup.find_max_overlap_indices(
                driven_uncoupled_state, return_args_coup['modes'][idx], return_overlaps=True
            )
            overlaps_coup_with_uncoupled.append(overlaps[0])
            energies_coup_with_uncoupled.append(return_args_coup['energies'][idx][overlap_idxs[0]])

            # Find the states in driven coupled modes closest to the dressed states
            overlap_idxs, overlaps = floquet_class_coup.find_max_overlap_indices(
                coup_ham.dressed_states[coup_ham.basis_to_dressed_mapping_dict[coup_ham.ket_to_basis_index(state_label)]],
                return_args_coup['modes'][idx],
                return_overlaps=True
            )
            overlaps_coup_with_dressed.append(overlaps[0])
            energies_coup_with_dressed.append(return_args_coup['energies'][idx][overlap_idxs[0]])
            results = {
                'overlaps_coup_with_uncoupled': overlaps_coup_with_uncoupled,
                'energies_coup_with_uncoupled': energies_coup_with_uncoupled,
                'overlaps_coup_with_dressed': overlaps_coup_with_dressed,
                'energies_coup_with_dressed': energies_coup_with_dressed
            }
        return results
    
    def plot_coupled_uncoupled_dressed_overlaps_and_energies(self,return_args_coup, 
                                                                 alice_results, 
                                                                 bob_results, 
                                                                 freq_unit=2 * np.pi, 
                                                                 energy_unit=2*np.pi):
        """
        Plot energies and overlaps for Alice and Bob with uncoupled and dressed states.

        Parameters
        ----------
        return_args_coup : dict
            Dictionary containing 'frequencies' (array of drive frequencies).
        alice_results : dict
            Dictionary with keys 'overlaps_coup_with_uncoupled', 'energies_coup_with_uncoupled',
            'overlaps_coup_with_dressed', 'energies_coup_with_dressed' for Alice.
        bob_results : dict
            Same as alice_results, but for Bob.
        freq_unit : float, optional
            Unit to divide frequencies by for plotting (default: 1e6 for MHz).
        energy_unit : float, optional
            Unit to divide energies by for plotting (default: 1e6 for MHz).
        """
        import matplotlib.pyplot as plt

        freqs = np.array(return_args_coup['frequencies']) / freq_unit

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        # Plot energies
        ax1.plot(freqs, np.array(alice_results['energies_coup_with_uncoupled']) / energy_unit, 'o-', label='Alice Energies with Uncoupled States')
        ax1.plot(freqs, np.array(alice_results['energies_coup_with_dressed']) / energy_unit, 'o-', label='Alice Energies with Dressed States')
        ax1.plot(freqs, np.array(bob_results['energies_coup_with_uncoupled']) / energy_unit, 'o-', label='Bob Energies with Uncoupled States')
        ax1.plot(freqs, np.array(bob_results['energies_coup_with_dressed']) / energy_unit, 'o-', label='Bob Energies with Dressed States')
        ax1.set_xlabel('Drive Frequency (MHz)')
        ax1.set_ylabel('Energy (MHz)')
        ax1.legend()
        ax1.grid()

        # Plot overlaps
        ax2.plot(freqs, alice_results['overlaps_coup_with_uncoupled'], 'o-', label='Alice Overlaps with Uncoupled States')
        ax2.plot(freqs, alice_results['overlaps_coup_with_dressed'], 'o-', label='Alice Overlaps with Dressed States')
        ax2.plot(freqs, bob_results['overlaps_coup_with_uncoupled'], 'o-', label='Bob Overlaps with Uncoupled States')
        ax2.plot(freqs, bob_results['overlaps_coup_with_dressed'], 'o-', label='Bob Overlaps with Dressed States')
        # draw a vertical line where alice uncoupled has minima
        min_idx = np.argmin(alice_results['overlaps_coup_with_uncoupled'])
        min_freq = freqs[min_idx]
        ax2.axvline(min_freq, color='k', linestyle='--', label='Alice Min Overlap')
        ax2.set_xlabel('Drive Frequency (MHz)')
        ax2.set_ylabel('Overlap')
        ax2.legend()
        ax2.grid()

        plt.tight_layout()
        plt.show()
        return min_freq, min_idx
    
    def calculate_overlaps_and_energies_plus_minus(self, floquet_class_uncoup, floquet_class_coup, return_args_uncoup, return_args_coup, state_label_alice, state_label_bob):
        """
        Calculates the overlaps and energies of the symmetric (plus) and antisymmetric (minus) combinations 
        of two quantum states (Alice and Bob) in both uncoupled and coupled Floquet systems.
        For each frequency in the coupled system, this method:
            - Identifies the basis states for Alice and Bob.
            - Finds the closest matching states in the uncoupled Floquet modes.
            - Constructs symmetric (plus) and antisymmetric (minus) superpositions of these states.
            - Computes overlaps and energies for these superpositions with respect to both the uncoupled and dressed (bare) states in the coupled system.
        Args:
            floquet_class_uncoup: An instance of the Floquet class for the uncoupled system, providing methods to find overlaps.
            floquet_class_coup: An instance of the Floquet class for the coupled system, providing methods to find overlaps.
            return_args_uncoup (dict): Dictionary containing results from the uncoupled Floquet calculation, 
                including 'modes' (list of mode arrays) and other relevant data.
            return_args_coup (dict): Dictionary containing results from the coupled Floquet calculation, 
                including 'modes' (list of mode arrays), 'energies' (list of energy arrays), and 'frequencies'.
            state_label_alice: Label or index identifying Alice's state in the basis.
            state_label_bob: Label or index identifying Bob's state in the basis.
        Returns:
            dict: A dictionary with the following keys, each containing a list (one entry per frequency):
                - 'plus_overlaps_coup_with_uncoupled': Overlaps of the coupled plus state with the uncoupled plus state.
                - 'plus_energies_coup_with_uncoupled': Energies of the coupled plus state closest to the uncoupled plus state.
                - 'plus_overlaps_coup_with_dressed': Overlaps of the coupled plus state with the dressed plus state.
                - 'plus_energies_coup_with_dressed': Energies of the coupled plus state closest to the dressed plus state.
                - 'minus_overlaps_coup_with_uncoupled': Overlaps of the coupled minus state with the uncoupled minus state.
                - 'minus_energies_coup_with_uncoupled': Energies of the coupled minus state closest to the uncoupled minus state.
                - 'minus_overlaps_coup_with_dressed': Overlaps of the coupled minus state with the dressed minus state.
                - 'minus_energies_coup_with_dressed': Energies of the coupled minus state closest to the dressed minus state.
        """
        
        coup_ham = self
        plus_overlaps_coup_with_uncoupled = []
        plus_energies_coup_with_uncoupled = []
        plus_overlaps_coup_with_dressed = []
        plus_energies_coup_with_dressed = []

        minus_overlaps_coup_with_uncoupled = []
        minus_energies_coup_with_uncoupled = []
        minus_overlaps_coup_with_dressed = []
        minus_energies_coup_with_dressed = []



        for idx, freq in enumerate(return_args_coup['frequencies']):
            alice_state = coup_ham.basis_states[coup_ham.ket_to_basis_index(state_label_alice)]
            bob_state = coup_ham.basis_states[coup_ham.ket_to_basis_index(state_label_bob)]

            # Find the state in the uncoupled modes
            alice_overlap_idxs, alice_overlaps = floquet_class_uncoup.find_max_overlap_indices(
                alice_state, return_args_uncoup['modes'][idx], return_overlaps=True
            )
            alice_driven_uncoupled_state = return_args_uncoup['modes'][idx][alice_overlap_idxs[0]]

            bob_overlap_idxs, bob_overlaps = floquet_class_uncoup.find_max_overlap_indices(
                bob_state, return_args_uncoup['modes'][idx], return_overlaps=True
            )
            bob_driven_uncoupled_state = return_args_uncoup['modes'][idx][bob_overlap_idxs[0]]

            

            # All states 
            plus_driven_uncoupled_state = (alice_driven_uncoupled_state + bob_driven_uncoupled_state) / np.sqrt(2)
            minus_driven_uncoupled_state = (alice_driven_uncoupled_state - bob_driven_uncoupled_state) / np.sqrt(2)
            plus_dressed_state = (alice_state + bob_state) / np.sqrt(2)
            minus_dressed_state = (alice_state - bob_state) / np.sqrt(2)
            
            # PLUS 
            # Find the states in driven coupled modes closest to the uncoupled plus state
            plus_overlap_idxs_coup, plus_overlaps_coup = floquet_class_coup.find_max_overlap_indices(
                plus_driven_uncoupled_state, return_args_coup['modes'][idx], return_overlaps=True
            )
            plus_overlaps_coup_with_uncoupled.append(plus_overlaps_coup[0])
            plus_energies_coup_with_uncoupled.append(return_args_coup['energies'][idx][plus_overlap_idxs_coup[0]])

            # Find the states in driven coupled modes closest to the dressed plus state
            
            plus_overlap_idxs_dressed, plus_overlaps_dressed = floquet_class_coup.find_max_overlap_indices(
                plus_dressed_state, return_args_coup['modes'][idx], return_overlaps=True
            )
            plus_overlaps_coup_with_dressed.append(plus_overlaps_dressed[0])
            plus_energies_coup_with_dressed.append(return_args_coup['energies'][idx][plus_overlap_idxs_dressed[0]])

            # Minus state

            # Find the states in driven coupled modes closest to the uncoupled minus state
            minus_overlap_idxs_coup, minus_overlaps_coup = floquet_class_coup.find_max_overlap_indices(
                minus_driven_uncoupled_state, return_args_coup['modes'][idx], return_overlaps=True
            )
                
            minus_overlaps_coup_with_uncoupled.append(minus_overlaps_coup[0])
            minus_energies_coup_with_uncoupled.append(return_args_coup['energies'][idx][minus_overlap_idxs_coup[0]])

            # Find the states in driven coupled modes closest to the dressed minus state
            minus_overlap_idxs_dressed, minus_overlaps_dressed = floquet_class_coup.find_max_overlap_indices(
                minus_dressed_state, return_args_coup['modes'][idx], return_overlaps=True
            )
            minus_overlaps_coup_with_dressed.append(minus_overlaps_dressed[0])
            minus_energies_coup_with_dressed.append(return_args_coup['energies'][idx][minus_overlap_idxs_dressed[0]])
            

        return_dict = {
            'plus_overlaps_coup_with_uncoupled': plus_overlaps_coup_with_uncoupled,
            'plus_energies_coup_with_uncoupled': plus_energies_coup_with_uncoupled,
            'plus_overlaps_coup_with_dressed': plus_overlaps_coup_with_dressed,
            'plus_energies_coup_with_dressed': plus_energies_coup_with_dressed,
            'minus_overlaps_coup_with_uncoupled': minus_overlaps_coup_with_uncoupled,
            'minus_energies_coup_with_uncoupled': minus_energies_coup_with_uncoupled,
            'minus_overlaps_coup_with_dressed': minus_overlaps_coup_with_dressed,
            'minus_energies_coup_with_dressed': minus_energies_coup_with_dressed
        }
        return return_dict
    
    def plot_plus_minus_overlaps_and_energies(self, return_args_coup, return_args_uncoup, plus_minus_results, min_idx):
        """
        Plot energies and overlaps for plus and minus states with uncoupled and dressed states.

        Parameters
        ----------
        return_args_coup : dict
            Dictionary containing 'frequencies' (array of drive frequencies) for the coupled system.
        return_args_uncoup : dict
            Dictionary containing 'frequencies' (array of drive frequencies) for the uncoupled system.
        plus_minus_results : dict
            Dictionary with keys:
                'plus_energies_coup_with_uncoupled', 'plus_energies_coup_with_dressed',
                'minus_energies_coup_with_uncoupled', 'minus_energies_coup_with_dressed',
                'plus_overlaps_coup_with_uncoupled', 'plus_overlaps_coup_with_dressed',
                'minus_overlaps_coup_with_uncoupled', 'minus_overlaps_coup_with_dressed'
        min_idx : int
            Index of the minimum overlap (for annotation).
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        # Plot energies
        ax1.plot(return_args_coup['frequencies'] / 2 / np.pi, np.array(plus_minus_results['plus_energies_coup_with_uncoupled']) / 2 / np.pi, 'o-', label='Plus Energies with Uncoupled States')
        ax1.plot(return_args_coup['frequencies'] / 2 / np.pi, np.array(plus_minus_results['plus_energies_coup_with_dressed']) / 2 / np.pi, 'o-', label='Plus Energies with Dressed States')
        ax1.plot(return_args_coup['frequencies'] / 2 / np.pi, np.array(plus_minus_results['minus_energies_coup_with_uncoupled']) / 2 / np.pi, 'o-', label='Minus Energies with Uncoupled States')
        ax1.plot(return_args_coup['frequencies'] / 2 / np.pi, np.array(plus_minus_results['minus_energies_coup_with_dressed']) / 2 / np.pi, 'o-', label='Minus Energies with Dressed States')

        # In the uncoupled case, plot the energy difference between plus and minus states
        differences = np.array(plus_minus_results['plus_energies_coup_with_uncoupled']) / 2 / np.pi - np.array(plus_minus_results['minus_energies_coup_with_uncoupled']) / 2 / np.pi
        ax1.plot(return_args_uncoup['frequencies'] / 2 / np.pi, differences, 'o--', label='Energy Difference (Uncoupled States)')
        # Draw a vertical line at min_idx and annotate the energy difference at that point
        min_freq = return_args_coup['frequencies'][min_idx] / 2 / np.pi
        energy_diff_at_min = differences[min_idx]
        ax1.axvline(min_freq, color='k', linestyle='--', label=f'Min Diff: {energy_diff_at_min:.2f} MHz')
        ax1.set_xlabel('Drive Frequency (MHz)')
        ax1.set_ylabel('Energy (MHz)')
        ax1.legend()
        ax1.grid()

        # Plot overlaps
        ax2.plot(return_args_coup['frequencies'] / 2 / np.pi, plus_minus_results['plus_overlaps_coup_with_uncoupled'], 'o-', label='Plus Overlaps with Uncoupled States')
        ax2.plot(return_args_coup['frequencies'] / 2 / np.pi, plus_minus_results['plus_overlaps_coup_with_dressed'], 'o-', label='Plus Overlaps with Dressed States')
        ax2.plot(return_args_coup['frequencies'] / 2 / np.pi, plus_minus_results['minus_overlaps_coup_with_uncoupled'], 'o-', label='Minus Overlaps with Uncoupled States')
        ax2.plot(return_args_coup['frequencies'] / 2 / np.pi, plus_minus_results['minus_overlaps_coup_with_dressed'], 'o-', label='Minus Overlaps with Dressed States')
        ax2.set_xlabel('Drive Frequency (MHz)')
        ax2.set_ylabel('Overlap')
        ax2.legend()
        ax2.grid()

        plt.tight_layout()
        plt.show()
        return energy_diff_at_min
    # ---------------------End of finding effective g between states---------------------


    
    