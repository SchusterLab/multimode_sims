import numpy as np
import h5py
import os
import matplotlib.pyplot as plt
import qutip as qt
from scripts.sims_base import sims_base
from tqdm import tqdm

class FloquetHamiltonian(sims_base):
    def __init__(self, H, T=None, args=None,  save_dir=None):
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
        super().__init__(sim_name = 'FloquetHamiltonian', save_dir=save_dir)
        self.H = H
        self.T = T
        if T is None:
            T = 2 * np.pi / args['omega']
        self.args = args
        self.options = qt.Options(nsteps=10000)
        # self.floquet_basis = qt.FloquetBasis(H, T, args, options=self.options)
    
    def init_floquet_basis(self):
        """
        Initialize the Floquet basis for the Hamiltonian.

        This method should be called before computing Floquet modes.
        """
        if not hasattr(self, 'H'):
            raise AttributeError("Hamiltonian (H) must be defined before initializing Floquet basis.")
        if not hasattr(self, 'T'):
            raise AttributeError("Floquet period (T) must be defined before initializing Floquet basis.")
        if not hasattr(self, 'args'):
            raise AttributeError("Arguments (args) must be defined before initializing Floquet basis.")
        
        self.floquet_basis = qt.FloquetBasis(self.H, self.T, self.args, options=self.options)

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

        filename = self.save_results_to_hdf5({'energies': quasienergies, 'modes': floquet_modes})

        return filename
    
    def load_floquet_modes(self, filename):
        """
        Load Floquet modes and quasienergies from an HDF5 file.

        Parameters
        ----------
        filename : str
            Name of the HDF5 file to load.

        Returns
        -------
        quasienergies : np.ndarray
            Array of quasienergies.
        floquet_modes : list of Qobj
            List of Floquet mode states.
        """
        file_path = os.path.join(self.save_dir, filename)
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"HDF5 file not found: {file_path}")

        with h5py.File(file_path, "r") as f:
            quasienergies = f['energies'][:]
            floquet_modes = [qt.Qobj(mode) for mode in f['modes'][:]]

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
        f_names = []

        for args in tqdm(args_list, desc="Computing Floquet modes"):
            self.args = args  # Update args for each computation
            T = 2 * np.pi / args['omega'] if 'omega' in args else self.T
            self.floquet_basis = qt.FloquetBasis(self.H, T, args, options=self.options)
            f_name = self.compute_floquet_modes()
            f_names.append(f_name)
            

        results = {
            # 'args_list': args_list,
            'filenames': f_names,
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
        fname = self.save_results_to_csv(results, filename='sweep_ac_drive_amplitude.csv')
        return f_name

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
        fname = self.save_results_to_csv(results, filename='sweep_ac_drive_frequency')
        return fname
    
    def load_sweep_results(self, filename, swept_param='frequencies'):
        '''
        use load_results_from_csv(f_name)
        '''
        fnames_dict = self.load_results_from_csv(filename)
        filenames = fnames_dict['filenames']
        swept_list = np.array(fnames_dict[swept_param])

        modes_swept_list = []
        energies_swept_list = []
        for i, f_name in enumerate(filenames):
            energies, modes = self.load_floquet_modes(f_name)
            modes_swept_list.append(modes)
            energies_swept_list.append(energies)
        results = {
            swept_param: swept_list,
            'modes': modes_swept_list,
            'energies': energies_swept_list}
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
           

    # mapping quasienergy to dressed energy 
    # fd approx n*omegad + fq ; figure out n that makes n*omega + fq closest to the dressed energy 
    # given fd and fq  
    def map_quasienergy_to_dressed_energy(self, quasienergy, dressed_energy, omega, tol=200* 2* np.pi):
        """
        Map a single dressed energy to the closest Floquet quasienergy shifted by n*omega.

        Parameters
        ----------
        quasienergy : float
            Floquet quasienergy.
        dressed_energy : float
            Dressed energy to map.
        omega : float
            Floquet drive frequency.
        tol : float, optional
            Tolerance for matching energies, in MHz. Default is 50.

        Returns
        -------
        result : dict
            Dictionary containing:
                - 'dressed_energy': The dressed energy.
                - 'closest_quasienergy': The closest quasienergy (shifted).
                - 'n': Integer n such that quasienergy + n*omega ≈ dressed_energy.
                - 'delta': Difference between dressed_energy and (quasienergy + n*omega).
                - 'matched': True if match within tolerance, else False.
        """
        n = int(np.round((dressed_energy - quasienergy) / omega))
        fq_shifted = quasienergy + n * omega
        delta = np.abs(dressed_energy - fq_shifted)
        matched = delta < tol
        return {
            'dressed_energy': dressed_energy,
            'closest_quasienergy': fq_shifted if matched else None,
            'n': n if matched else None,
            'delta': delta,
            'matched': matched
        }
