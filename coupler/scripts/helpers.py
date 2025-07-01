import numpy as np
import qutip as qt
from tqdm import tqdm
import matplotlib.pyplot as plt
from scripts.sims_base import sims_base
from scripts.floquet import *
import sys
from IPython.display import clear_output
import math

# There should be 2 classes , one that is just solving floquet/hamiltonian diagnalization. Let 
# call it FLoquetHamiltonian class. And the second class is its daughter class, which solves it 
# for our particular problem that involves 3 objects and couper



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
        self.generate_basis_energies(self.generate_H0(coupling=False))
        self.generate_dressed_states()
    
    def get_participation_factors(self): 
        """
        Calculate the participation factors for the dressed states.

        Returns
        -------
        participation_factors : list of float
            List of participation factors for each dressed state.
        """
        participation_factors = {}
        alice_delta = self.omega_A - self.omega_C
        bob_delta = self.omega_B - self.omega_C
        alice_part = self.g_AC / alice_delta if alice_delta != 0 else 0
        bob_part = self.g_BC / bob_delta if bob_delta != 0 else 0
        participation_factors = {
            'alice': alice_part,
            'bob': bob_part,}
        return participation_factors
    
    def get_dressed_energy(self, ket_label):
        """
        Get the dressed energy for a given ket label.

        Parameters
        ----------
        ket_label : str
            Ket label in the form '|n_A,n_B,n_C>'.

        Returns
        -------
        energy : float
            The dressed energy corresponding to the ket label.
        """
        index = self.basis_to_dressed_mapping_dict[self.ket_to_basis_index(ket_label)]
        print(f"Index for ket {ket_label}: {index}")
        return self.dressed_energies[index]

    def get_dressed_state(self, ket_label):
        """
        Get the dressed state for a given ket label.

        Parameters
        ----------
        ket_label : str
            Ket label in the form '|n_A,n_B,n_C>'.

        Returns
        -------
        dressed_state : Qobj
            The dressed state corresponding to the ket label.
        """
        index = self.basis_to_dressed_mapping_dict[self.ket_to_basis_index(ket_label)]
        return self.dressed_states[index]
    def get_basis_state(self, ket_label):
        """
        Get the basis state for a given ket label.

        Parameters
        ----------
        ket_label : str
            Ket label in the form '|n_A,n_B,n_C>'.

        Returns
        -------
        basis_state : Qobj
            The basis state corresponding to the ket label.
        """
        index = self.ket_to_basis_index(ket_label)
        return self.basis_states[index]

    def generate_basis_operators(self):
        I = qt.qeye(self.trunc)
        a = qt.tensor(qt.destroy(self.trunc), I, I)
        a_dag = a.dag()
        b = qt.tensor(I, qt.destroy(self.trunc), I)
        b_dag = b.dag()
        c = qt.tensor(I, I, qt.destroy(self.trunc))
        c_dag = c.dag()
        return a, a_dag, b, b_dag, c, c_dag
    


    def generate_H0(self, coupling = True):
        # if rotating:
        #     return qt.qeye(self.trunc ** 3)
        
        H = (
            self.omega_A * self.a_dag * self.a
            + self.omega_B * self.b_dag * self.b
            + self.omega_C * self.c_dag * self.c
        )
        if not coupling:
            print("No couplings added to Hamiltonian")
        else:
            print("Adding jc couplings to Hamiltonian")
            H += (
                # self.g_AC * (self.a_dag *self.c + self.a * self.c_dag)#
                self.g_AC *(self.a_dag - self.a) * (self.c_dag - self.c)
               +  self.g_BC * (self.b_dag - self.b) * (self.c_dag - self.c)
            )
        
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
        self.theta_zpt = theta_zpt
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
    
    # Generate static modes (no coupling)
    def generate_basis_energies(self, H_static):
        """
        Diagonalize the static Hamiltonian H_static and return the eigenvalues and eigenstates.

        Parameters
        ----------
        H_static : Qobj
            The static Hamiltonian to diagonalize.

        Returns
        -------
        eigenvalues : np.ndarray
            Array of eigenvalues.
        eigenstates : list of Qobj
            List of eigenstate Qobjs.
        """
        eigenvalues, eigenstates = H_static.eigenstates()
        # eigenstates contain same states as self.basis_states but in different order
        # can you reorder eigenvalues to match the order of self.basis_states?
        # Reorder eigenvalues to match the order of self.basis_states
        # For each basis state, find which eigenstate it overlaps with most, and reorder accordingly
        reordered_eigenvalues = np.zeros_like(eigenvalues)
        reordered_eigenstates = [None] * len(eigenstates)
        for i, basis_state in enumerate(self.basis_states):
            overlaps = [np.abs(basis_state.overlap(eig_state))**2 for eig_state in eigenstates]
            max_idx = int(np.argmax(overlaps))
            reordered_eigenvalues[i] = eigenvalues[max_idx]
            reordered_eigenstates[i] = eigenstates[max_idx]
        eigenvalues = reordered_eigenvalues
        eigenstates = reordered_eigenstates

        self.basis_energies = eigenvalues
        return None
    
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
        self.H0 = self.generate_H0( coupling=True)
        operator, time_func = self.linc_potential_operator()
        H_static = self.H0 #+ (time_func(0, args={'amp': 0, 'omega': 0}) * operator)

        dressed_static_energies, dressed_static_states = self.static_modes(H_static)
        basis_to_dressed_mapping_dict, dressed_to_basis_mapping_dict = self.basis_to_dressed_mapping(dressed_static_states)

        self.dressed_energies = dressed_static_energies
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

    def plot_basis_state_overlap(self, dressed_state, max_total_photons=3, title=None, yscale = 'linear'):
        """
        Plot the overlap of a dressed state with the basis states, ignoring states with more than max_total_photons.

        Parameters
        ----------
        dressed_state : Qobj
            The dressed state to analyze.
        max_total_photons : int, optional
            Maximum total number of photons allowed in a basis state. Default is 3.
        title : str, optional
            Title for the plot.

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
        if title is not None:
            plt.title(title)
        else:
            plt.title('Overlap of Dressed State with Basis States (Filtered by Total Photons)')
        plt.xticks(filtered_indices, [self.basis_index_to_ket(i) for i in filtered_indices], rotation=90)
        plt.yscale(yscale)
        if yscale == 'log':
            plt.ylim(bottom=1e-3)
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
