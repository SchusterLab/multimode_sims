import numpy as np
import qutip as qt
from tqdm import tqdm
import matplotlib.pyplot as plt
from scripts.sims_base import sims_base
from scripts.floquet import *
import sys
from IPython.display import clear_output
import math

class coupler_driving_analysis(sims_base):
    def __init__(self, sim_name, save_dir, plots_dir, coupler_ham=None): 
        """
        Initialize the coupler_driving_analysis class.

        Parameters
        ----------
        sim_name : str
            Name of the simulation.
        dir_name : str
            Directory name for saving results.
        """
        super().__init__(sim_name, save_dir, plots_dir)
        self.coupler_ham = coupler_ham
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
        idxs_coup_with_uncoupled = []
        overlaps_coup_with_dressed = []

        energies_coup_with_dressed = []
        idxs_coup_with_dressed = []
        coup_ham = self.coupler_ham
        prev_idx = None 

        for idx, freq in enumerate(return_args_coup['frequencies']):
            state = coup_ham.basis_states[coup_ham.ket_to_basis_index(state_label)]

            # Find the state in the driven uncoupled modes
            overlap_idxs_uncoup, overlaps_uncoup = floquet_class_uncoup.find_max_overlap_indices(
            state, return_args_uncoup['modes'][idx], return_overlaps=True
            )
            driven_uncoupled_state = return_args_uncoup['modes'][idx][overlap_idxs_uncoup[0]]

            # Find the driven coupled state with lagrest overlap with driven uncoupled state 
            overlap_idxs_coup_uncoup, overlaps_coup_uncoup = floquet_class_coup.find_max_overlap_indices(
            driven_uncoupled_state, return_args_coup['modes'][idx], return_overlaps=True
            )
            overlaps_coup_with_uncoupled.append(overlaps_coup_uncoup[0])
            energies_coup_with_uncoupled.append(return_args_coup['energies'][idx][overlap_idxs_coup_uncoup[0]])
            idxs_coup_with_uncoupled.append(overlap_idxs_coup_uncoup[0])

            # Find the states in driven coupled modes closest to the dressed state
            dressed_state = coup_ham.dressed_states[coup_ham.basis_to_dressed_mapping_dict[coup_ham.ket_to_basis_index(state_label)]]
            overlap_idxs_coup_dressed, overlaps_coup_dressed = floquet_class_coup.find_max_overlap_indices(
            dressed_state,
            return_args_coup['modes'][idx],
            return_overlaps=True
            )

            # when plotting energy as function of frequency, sometimes the identified energy does not fit (causes a kink in the graph)
            # a check to prevent this if there are two states which are very close in overlap, choose the one with the closest energy to the previous idx

            
            overlaps_coup_with_dressed.append(overlaps_coup_dressed[0])
            energies_coup_with_dressed.append(return_args_coup['energies'][idx][overlap_idxs_coup_dressed[0]])
            idxs_coup_with_dressed.append(overlap_idxs_coup_dressed[0])
            prev_idx = overlap_idxs_coup_dressed[0]

        results = {
            'overlaps_coup_with_uncoupled': overlaps_coup_with_uncoupled,
            'energies_coup_with_uncoupled': energies_coup_with_uncoupled,
            'idxs_coup_with_uncoupled': idxs_coup_with_uncoupled,
            'overlaps_coup_with_dressed': overlaps_coup_with_dressed,
            'energies_coup_with_dressed': energies_coup_with_dressed,
            'idxs_coup_with_dressed': idxs_coup_with_dressed
        }
        return results

    def calculate_coupled_uncoupled_dressed_overlaps_and_energies_two_photon(self, floquet_class_uncoup, floquet_class_coup,  return_args_uncoup, return_args_coup, 
         two_photon_state_label, single_photon_results, dag_op):           
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

        06/17/2025: for 2 photon state, apply adag to corresponding single photon floquet state, find overlaps with this 
        adaged state .
        """
        
        overlaps_coup_with_uncoupled = []
        energies_coup_with_uncoupled = []
        idxs_coup_with_uncoupled = []
        overlaps_coup_with_dressed = []

        energies_coup_with_dressed = []
        idxs_coup_with_dressed = []
        coup_ham = self.coupler_ham
        prev_idx = None 

        # print(f"Calculating overlaps and energies for two-photon state: {two_photon_state_label}")
        # print(single_photon_results['idxs_coup_with_uncoupled'])

        for idx, freq in enumerate(return_args_coup['frequencies']):
            # print(f"Processing frequency index {idx} with frequency {freq / (2 * np.pi)} MHz")
            state = coup_ham.basis_states[coup_ham.ket_to_basis_index(two_photon_state_label)]

            
            # Find the state in the driven uncoupled adaged state 
            # print(single_photon_results['idxs_coup_with_uncoupled'][idx])
            single_photon_floquet_state = return_args_coup['modes'][idx][single_photon_results['idxs_coup_with_uncoupled'][idx]]
            single_photon_floquet_state_formatted =qt.Qobj(single_photon_floquet_state.full().flatten(), dims=state.dims) # format it 
            adaged_driven_uncoupled_state = (dag_op * single_photon_floquet_state_formatted).unit()
            

            # Find the driven coupled state with lagrest overlap with driven uncoupled state 
            overlap_idxs_coup_uncoup, overlaps_coup_uncoup = floquet_class_coup.find_max_overlap_indices(
            adaged_driven_uncoupled_state, return_args_coup['modes'][idx], return_overlaps=True
            )
            overlaps_coup_with_uncoupled.append(overlaps_coup_uncoup[0])
            energies_coup_with_uncoupled.append(return_args_coup['energies'][idx][overlap_idxs_coup_uncoup[0]])
            idxs_coup_with_uncoupled.append(overlap_idxs_coup_uncoup[0])

            # Find the states in driven coupled modes closest to the dressed state
            # Find the state in the dressed ststae  
            # print('now dressed')
            single_photon_floquet_state = return_args_coup['modes'][idx][single_photon_results['idxs_coup_with_dressed'][idx]]
            single_photon_floquet_state_formatted =qt.Qobj(single_photon_floquet_state.full().flatten(), dims=state.dims) # format it 
            adaged_floquet_state = (dag_op * single_photon_floquet_state_formatted).unit()
            
            overlap_idxs_coup_dressed, overlaps_coup_dressed = floquet_class_coup.find_max_overlap_indices(
            adaged_floquet_state,
            return_args_coup['modes'][idx],
            return_overlaps=True
            )

            # when plotting energy as function of frequency, sometimes the identified energy does not fit (causes a kink in the graph)
            # a check to prevent this if there are two states which are very close in overlap, choose the one with the closest energy to the previous idx

            
            overlaps_coup_with_dressed.append(overlaps_coup_dressed[0])
            energies_coup_with_dressed.append(return_args_coup['energies'][idx][overlap_idxs_coup_dressed[0]])
            idxs_coup_with_dressed.append(overlap_idxs_coup_dressed[0])
            prev_idx = overlap_idxs_coup_dressed[0]

        results = {
            'overlaps_coup_with_uncoupled': overlaps_coup_with_uncoupled,
            'energies_coup_with_uncoupled': energies_coup_with_uncoupled,
            'idxs_coup_with_uncoupled': idxs_coup_with_uncoupled,
            'overlaps_coup_with_dressed': overlaps_coup_with_dressed,
            'energies_coup_with_dressed': energies_coup_with_dressed,
            'idxs_coup_with_dressed': idxs_coup_with_dressed
        }
        return results
    
    def plot_coupled_uncoupled_dressed_overlaps_and_energies(self, return_args_coup, 
                                                                 alice_results, 
                                                                 bob_results, 
                                                                 freq_unit=2 * np.pi, 
                                                                 energy_unit=2*np.pi,
                                                                 save=True,
                                                                 save_fig=True):
        """
        Plot energies, overlaps, and dressed overlap indices for Alice and Bob.

        Parameters
        ----------
        return_args_coup : dict
            Dictionary containing 'frequencies' (array of drive frequencies).
        alice_results : dict
            Dictionary with keys 'overlaps_coup_with_uncoupled', 'energies_coup_with_uncoupled',
            'overlaps_coup_with_dressed', 'energies_coup_with_dressed', 'idxs_coup_with_dressed' for Alice.
        bob_results : dict
            Same as alice_results, but for Bob.
        freq_unit : float, optional
            Unit to divide frequencies by for plotting (default: 1e6 for MHz).
        energy_unit : float, optional
            Unit to divide energies by for plotting (default: 1e6 for MHz).
        save_fig : bool, optional
            If True, saves the figure.
        """
        import matplotlib.pyplot as plt

        freqs = np.array(return_args_coup['frequencies']) / freq_unit

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))

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

        # Plot dressed overlap indices
        ax3.plot(freqs, alice_results.get('idxs_coup_with_dressed', []), 'o-', label='Alice Dressed Overlap Idx')
        ax3.plot(freqs, bob_results.get('idxs_coup_with_dressed', []), 'o-', label='Bob Dressed Overlap Idx')
        ax3.axvline(min_freq, color='k', linestyle='--', label='Alice Min Overlap')
        ax3.set_xlabel('Drive Frequency (MHz)')
        ax3.set_ylabel('Dressed Overlap Index')
        ax3.legend()
        ax3.grid()

        plt.tight_layout()
        plt.show()
        if save or save_fig:
            self.save_plot_and_log(fig, 'coupled_uncoupled_dressed_overlaps_and_energies')
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
        
        coup_ham = self.coupler_ham
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
    
    def find_local_minimum_in_window(self, differences, frequencies, center_idx, window):
        """
        Find the local minimum of 'differences' within a window around 'center_idx'.

        Parameters
        ----------
        differences : array-like
            Array of values (e.g., energy differences).
        frequencies : array-like
            Array of frequencies corresponding to the differences.
        center_idx : int
            Center index around which to search for the local minimum.
        window : int
            Number of indices on each side of center_idx to include in the window.

        Returns
        -------
        local_min_idx : int
            Index of the local minimum within the window (relative to the full array).
        min_freq : float
            Frequency at the local minimum.
        min_diff : float
            Value of 'differences' at the local minimum.
        """
        start = max(center_idx - window, 0)
        end = min(center_idx + window + 1, len(differences))
        window_diffs = differences[start:end]
        local_min_idx_in_window = np.argmin(np.abs(window_diffs))
        local_min_idx = start + local_min_idx_in_window
        min_freq = frequencies[local_min_idx] / 2 / np.pi
        min_diff = differences[local_min_idx]
        return local_min_idx, min_freq, min_diff

    def plot_plus_minus_overlaps_and_energies(self, return_args_coup, return_args_uncoup, plus_minus_results, min_idx, return_all_diffs = False,
                                              save = True, save_fig=True):
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
        save_fig : bool, optional
            If True, saves the figure.
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
        # Find the minimum energy difference in a window of ±5 indices around min_idx
        window = 10
        local_min_idx, min_freq, local_energy_diff_at_min = self.find_local_minimum_in_window(
            differences, return_args_coup['frequencies'], min_idx, window
        )
        ax1.axvline(min_freq, color='r', linestyle='--', label=f'Local Min Diff: {local_energy_diff_at_min:.2f} MHz')
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
        if save or save_fig:
            self.save_plot_and_log(fig, 'coupled_uncoupled_dressed_overlaps_and_energies')
        
        return local_energy_diff_at_min
    
    def sweep_floquet_amplitude_and_analyze(self, amp_list, min_freq_for_fq, max_freq_for_fq, freq_pts, args=None, plot=True
        ):
        """
            Sweeps over a list of Floquet drive amplitudes, analyzes the coupled and uncoupled system dynamics,
            and extracts relevant physical quantities such as minimum overlap frequencies and energy differences.
            For each amplitude in `amp_list`, this method:
                - Constructs the uncoupled and coupled Hamiltonians using `coup_ham`.
                - Performs a frequency sweep for each case using the Floquet formalism.
                - Calculates overlaps and energies for specific basis states (e.g., '|1,0,0>', '|0,1,0>').
                - Optionally plots the results or finds minima analytically.
                - Computes and stores the minimum frequency (delta_bs) and minimum energy difference (g_bs).
                - Collects detailed results for further analysis.
            Args:
                coup_ham: An object providing methods to generate Hamiltonians and perform analysis.
                amp_list (list or array): List of Floquet drive amplitudes to sweep over.
                min_freq_for_fq (float): Minimum frequency for the Floquet frequency sweep.
                max_freq_for_fq (float): Maximum frequency for the Floquet frequency sweep.
                freq_pts (int): Number of frequency points in the sweep.
                args (dict): Additional arguments passed to the FloquetHamiltonian.
                plot (bool, optional): If True, plots analysis results; otherwise, only computes and returns values. Default is True.
            Returns:
                tuple:
                    delta_bs (list): List of minimum overlap frequencies for each amplitude.
                    g_bs (list): List of minimum energy differences for each amplitude.
                    results (list): List of dictionaries containing detailed analysis results for each amplitude.
        """

        delta_bs = []
        g_bs = []
        coup_ham = self.coupler_ham
        uncoup_fnames = []
        coup_fnames = []

        for idx, amp in enumerate(amp_list):
            

            # uncoupled case 
            H0 = coup_ham.generate_H0(coupling=False)
            operator, time_func = coup_ham.linc_potential_operator()
            H_fq = [H0, [operator, time_func]]
            floquet_class_uncoup = FloquetHamiltonian(H=H_fq, T=2*np.pi/200, args=args, save_dir=self.save_dir)
            f_name = floquet_class_uncoup.sweep_ac_drive_frequency(
                amp, freq_range=(min_freq_for_fq, max_freq_for_fq), num_points=freq_pts
            )
            return_args_uncoup = floquet_class_uncoup.load_sweep_results(f_name, 'frequencies')
            uncoup_fnames.append(f_name)

            # coupled case
            H0 = coup_ham.generate_H0(coupling=True)
            operator, time_func = coup_ham.linc_potential_operator()
            H_fq = [H0, [operator, time_func]]
            floquet_class_coup = FloquetHamiltonian(H=H_fq, T=2*np.pi/200, args=args, save_dir=self.save_dir)
            f_name = floquet_class_coup.sweep_ac_drive_frequency(
                amp, freq_range=(min_freq_for_fq, max_freq_for_fq), num_points=freq_pts
            )
            return_args_coup = floquet_class_coup.load_sweep_results(f_name, 'frequencies')
            coup_fnames.append(f_name)

            # Analysis
            analysis_dict = self.run_full_g_kerr_analysis(
                floquet_class_uncoup, floquet_class_coup, return_args_uncoup, return_args_coup)
            g_bs.append(analysis_dict['g_BS'])
            delta_bs.append(analysis_dict['delta_BS'])
            # Close previous plots
            plt.close('all')
            clear_output(wait=True)
            self.plot_floquet_amp_sweep_results(
                amp_list[:len(delta_bs)], delta_bs, g_bs, save = False
            )
        
        results = {
            'amp_list': amp_list,
            'uncoupled_fnames': uncoup_fnames,
            'coupled_fnames': coup_fnames,
            'delta_bs': delta_bs,
            'g_bs': g_bs,

        }
        self.save_results_to_csv(results, filename='sweep_amps_and_frequencies.csv')
        self.plot_floquet_amp_sweep_results(amp_list, delta_bs, g_bs, save=True)

        return results  
    
    ### Do analysis 
    def run_full_g_kerr_analysis(self, floquet_class_uncoup, floquet_class_coup, return_args_uncoup, return_args_coup, save_fig = True):
        """
        Run full g and Kerr analysis: calculates overlaps, energies, frequency shifts, and Kerr shifts.

        Parameters
        ----------
        floquet_class_uncoup : object
            FloquetHamiltonian instance for the uncoupled system.
        floquet_class_coup : object
            FloquetHamiltonian instance for the coupled system.
        return_args_uncoup : dict
            Results from uncoupled Floquet sweep.
        return_args_coup : dict
            Results from coupled Floquet sweep.
        min_idx : int
            Index of minimum overlap (for annotation).

        Returns
        -------
        results : dict
            Dictionary containing all analysis results, organized by category.
        """
        # floquet class not important
        if floquet_class_uncoup is None:
            H0_uncoup = self.coupler_ham.generate_H0(no_coupling=True)
            operator, time_func = self.coupler_ham.linc_potential_operator()
            H_fq_uncoup = [H0_uncoup, [operator, time_func]]
            floquet_class_uncoup = FloquetHamiltonian(H=H_fq_uncoup, T=2*np.pi/200, args=None, save_dir=self.save_dir)

        if floquet_class_coup is None:
            H0_coup = self.coupler_ham.generate_H0(no_coupling=False)
            operator, time_func = self.coupler_ham.linc_potential_operator()
            H_fq_coup = [H0_coup, [operator, time_func]]
            floquet_class_coup = FloquetHamiltonian(H=H_fq_coup, T=2*np.pi/200, args=None, save_dir=self.save_dir)

        print("Running full g and Kerr analysis...")
        # single_photon_floquet_state = return_args_coup['modes'][0][25]
        # print(f"single photon floquet state: {single_photon_floquet_state}")

        coup_ham = self.coupler_ham
        state_labels = ['|0,0,0>', '|1,0,0>', '|0,1,0>',  '|0,0,1>']
        results_dict = {}
        for state in state_labels:
            results_dict[state] = self.calculate_coupled_uncoupled_dressed_overlaps_and_energies(
                floquet_class_uncoup, floquet_class_coup, return_args_uncoup, return_args_coup, state
            )
        # for two photon states 
        two_photon_state_labels = ['|2,0,0>', '|0,2,0>', '|0,0,2>']
        corresponding_single_photon_results = {
            '|2,0,0>': results_dict['|1,0,0>'],
            '|0,2,0>': results_dict['|0,1,0>'],
            '|0,0,2>': results_dict['|0,0,1>']
        }
        corresponding_operators= {
            '|2,0,0>': coup_ham.a_dag,
            '|0,2,0>': coup_ham.b_dag,
            '|0,0,2>': coup_ham.c_dag
        }
        for two_photon_state_label in two_photon_state_labels: 
            # get the corresponding single photon results 
            single_photon_results = corresponding_single_photon_results[two_photon_state_label]
            adag_op = corresponding_operators[two_photon_state_label]
            results_dict[two_photon_state_label] = self.calculate_coupled_uncoupled_dressed_overlaps_and_energies_two_photon(
                floquet_class_uncoup, floquet_class_coup, return_args_uncoup, return_args_coup, two_photon_state_label,
                single_photon_results, adag_op
            )
            print(f"Processed two-photon state: {two_photon_state_label}")
        
        all_state_labels = state_labels + two_photon_state_labels


        # get the g and the bs     
        alice_results = results_dict['|1,0,0>']
        bob_results = results_dict['|0,1,0>']
        min_freq, min_idx = self.plot_coupled_uncoupled_dressed_overlaps_and_energies(
                return_args_coup, alice_results, bob_results, save = False, save_fig=save_fig
            )
        print(f"Minimum frequency: {min_freq} MHz, Minimum index: {min_idx}")


        plus_minus_results = self.calculate_overlaps_and_energies_plus_minus(
            floquet_class_uncoup, floquet_class_coup, return_args_uncoup, return_args_coup, '|1,0,0>', '|0,1,0>'
        )
        min_energy_diff = self.plot_plus_minus_overlaps_and_energies(
            return_args_coup, return_args_uncoup, plus_minus_results, min_idx, save=False, save_fig=save_fig
        )

        # get the frequency shift 
        freq_shifts, dressed_energies = self.plot_dressed_state_frequency_shift(
            all_state_labels, return_args_coup, floquet_class_coup, results_dict, save=save_fig
        )
        # get the kerr
        kerr_dict = self.plot_kerr_shifts(return_args_coup, freq_shifts, dressed_energies, save_fig=save_fig)

        results = {
            # 'overlap_energy_results': results_dict,
            # 'plus_minus_results': plus_minus_results,
            'g_BS': min_energy_diff /2,
            'delta_BS': min_freq,
            'freq_shifts': freq_shifts,
            'dressed_energies': dressed_energies,
            'kerr_dict': kerr_dict,
        }
        return results
    
    def get_kerr(self, no_photon_energy, single_photon_energy, two_photon_energy):
        """
        Calculate the Kerr shift from the energies of the states.

        Parameters
        ----------
        no_photon_energy : float
            Energy of the |0,0,0> state.
        singletone_energy : float
            Energy of the |1,0,0> state.
        two_photon_energy : float
            Energy of the |2,0,0> state.

        Returns
        -------
        kerr_shift : float
            The calculated Kerr shift.
        """
        return (two_photon_energy + no_photon_energy - 2 * single_photon_energy) 


    # Kerr shift function 
    def plot_kerr_shifts(self, return_args_coup, freq_shifts, dressed_energies, save_fig=True):
        """
        Plot the Kerr shifts for Alice, Bob, and Coupler as a function of drive frequency.

        Parameters
        ----------
        return_args_coup : dict
            Dictionary containing 'frequencies' for the coupled Floquet calculation.
        freq_shifts : dict
            Dictionary mapping state label to array of driven energies.
        dressed_energies : dict
            Dictionary mapping state label to dressed state energy (static, not driven).
        save_fig : bool, optional
            If True, saves the figure.

        Returns
        -------
        kerr_dict : dict
            Dictionary containing static and driven Kerr shifts for Alice, Bob, and Coupler.
        """
        alice_kerr = (freq_shifts['|2,0,0>'] + freq_shifts['|0,0,0>']) - 2 * freq_shifts['|1,0,0>']
        alice_static_kerr = (dressed_energies['|2,0,0>'] + dressed_energies['|0,0,0>']) - 2 * dressed_energies['|1,0,0>']
        bob_kerr = (freq_shifts['|0,2,0>'] + freq_shifts['|0,0,0>']) - 2 * freq_shifts['|0,1,0>']
        bob_static_kerr = (dressed_energies['|0,2,0>'] + dressed_energies['|0,0,0>']) - 2 * dressed_energies['|0,1,0>']
        coupler_kerr = (freq_shifts['|0,0,2>'] + freq_shifts['|0,0,0>']) - 2 * freq_shifts['|0,0,1>']
        coupler_static_kerr = (dressed_energies['|0,0,2>'] + dressed_energies['|0,0,0>']) - 2 * dressed_energies['|0,0,1>']

        plt.figure(figsize=(8, 6))
        plt.plot(return_args_coup['frequencies'] / 2 / np.pi, np.abs(alice_kerr / 2 / np.pi), label='Alice Kerr Shift', color='blue')
        plt.axhline(np.abs(alice_static_kerr / 2 / np.pi), color='blue', linestyle='--', label='Alice Static Kerr Shift')
        plt.plot(return_args_coup['frequencies'] / 2 / np.pi, np.abs(bob_kerr / 2 / np.pi), label='Bob Kerr Shift', color='orange')
        plt.axhline(np.abs(bob_static_kerr / 2 / np.pi), color='orange', linestyle='--', label='Bob Static Kerr Shift')
        plt.plot(return_args_coup['frequencies'] / 2 / np.pi, np.abs(coupler_kerr / 2 / np.pi), label='Coupler Kerr Shift', color='green')
        plt.axhline(np.abs(coupler_static_kerr / 2 / np.pi), color='green', linestyle='--', label='Coupler Static Kerr Shift')
        plt.xlabel('Drive Frequency (MHz)')
        plt.ylabel('Kerr Shift (MHz)')
        plt.title('Kerr Shift as a Function of Drive Frequency (log scale)')
        plt.yscale('log')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        if save_fig:
            self.save_plot_and_log(plt.gcf(), 'kerr_shifts_vs_drive_frequency')
        plt.show()

        kerr_dict = {
            'alice_kerr': alice_kerr,
            'alice_static_kerr': alice_static_kerr,
            'bob_kerr': bob_kerr,
            'bob_static_kerr': bob_static_kerr,
            'coupler_kerr': coupler_kerr,
            'coupler_static_kerr': coupler_static_kerr,
        }
        return kerr_dict
    


    def plot_floquet_amp_sweep_results(self, amp_list, delta_bs, g_bs, save = False, save_fig=True):
        fig, ax1 = plt.subplots(figsize=(8, 5))
        ax1.plot(amp_list, delta_bs, 'o-', label=r'$\Delta_{bs}$ (MHz)')
        ax1.set_xlabel('Drive Amplitude')
        ax2 = ax1.twinx()
        ax2.plot(amp_list, g_bs, 's-', color='orange', label=r'$g_{bs}$ (MHz)')
        ax1.set_xlabel(r'Drive Amplitude')
        ax1.set_ylabel(r'$\Delta_{bs}$ (MHz)')
        ax2.set_ylabel(r'$g_{bs}$ (MHz)')
        fig.legend(loc='upper left')
        plt.title('Sweep Results')
        plt.show()
        if save or save_fig:
            self.save_plot_and_log(fig, 'floquet_amp_sweep_results')
    

    ### Finding the frequency shift of a particular state when drive is on 
    # Given a state label, find the closest quasienergy to the dressed energy and plot it 
    # the function could be given multiple state lables, do for all 
    def plot_dressed_state_frequency_shift(self, state_labels, return_args_coup, floquet_class_coup, results_dict, plot=True, 
        save=True, save_fig=True):     
        """
        Plot the frequency shift of specified dressed states under drive.

        Parameters
        ----------
        state_labels : list of str
            List of ket strings for the states to analyze (e.g., ['|1,0,0>', '|0,1,0>']).
        return_args_coup : dict
            Dictionary containing 'frequencies' and 'energies' for the coupled Floquet calculation.
        floquet_class_coup : object
            FloquetHamiltonian instance for the coupled system, must have map_quasienergy_to_dressed_energy().
        results_dict : dict
            Dictionary with keys like 'energies_coup_with_uncoupled' for each state.
        plot : bool, optional
            If True, plot the results.
        save_fig : bool, optional
            If True, saves the figure.

        Returns
        -------
        freq_shifts : dict
            Dictionary mapping state label to array of driven energies.
        dressed_energies_dict : dict
            Dictionary mapping state label to dressed state energy (static, not driven).
        """

        freq_shifts = {}
        dressed_energies_dict = {}
        coup_ham = self.coupler_ham
        n_states = len(state_labels)
        ncols = min(3, n_states)
        nrows = math.ceil(n_states / ncols)

        if plot:
            fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False)

        for idx_state, state_label in enumerate(state_labels):
            dressed_idx = coup_ham.basis_to_dressed_mapping_dict[coup_ham.ket_to_basis_index(state_label)]
            dressed_energy = coup_ham.dressed_energies[dressed_idx]
            dressed_energies_dict[state_label] = dressed_energy
            driven_energies = []
            for idx, freq in enumerate(return_args_coup['frequencies']):
                quasienergy = results_dict[state_label]['energies_coup_with_uncoupled'][idx]
                mapping = floquet_class_coup.map_quasienergy_to_dressed_energy(
                    quasienergy, dressed_energy, freq
                )
                driven_energies.append(mapping['closest_quasienergy'])
            freq_shifts[state_label] = np.array(driven_energies)
            if plot:
                row = idx_state // ncols
                col = idx_state % ncols
                ax = axes[row][col]
                ax.plot(
                    return_args_coup['frequencies'] / 2 / np.pi,
                    freq_shifts[state_label] / 2 / np.pi,
                    label=f'{state_label} Driven Energy'
                )
                ax.axhline(dressed_energy / 2 / np.pi, color='gray', linestyle='--', label='Dressed Energy')
                ax.set_xlabel('Drive Frequency (MHz)')
                ax.set_ylabel('Driven State Energy (MHz)')
                ax.set_title(f'Driven State: {state_label}')
                ax.legend()
                ax.grid(True)

        if plot:
            # Hide unused subplots
            for idx in range(n_states, nrows * ncols):
                row = idx // ncols
                col = idx % ncols
                fig.delaxes(axes[row][col])
            plt.tight_layout()
            plt.show()
            if save or save_fig:
                self.save_plot_and_log(fig, 'dressed_state_frequency_shift')
        return freq_shifts, dressed_energies_dict


    
    