import numpy as np
import h5py
import os
import matplotlib.pyplot as plt
import qutip as qt
import pandas as pd
import datetime
import shutil # for moving files

# problem 06/10: some files stored in the save_dir, some in the plots_dir

class sims_base(): 
    def __init__(self, sim_name, save_dir, plots_dir=None):
        """
        Parameters
        ----------
        sim_name : str
            Name of the simulation.
        save_dir : str
            Directory where data/results will be saved.
        plots_dir : str, optional
            Directory where plots and logs will be saved. If None, defaults to save_dir.
        """
        self.sim_name = sim_name
        self.save_dir = save_dir
        self.plots_base_dir = plots_dir if plots_dir is not None else save_dir
        # Ensure logs for plots are written in the plots directory, not the save directory
        self.log_dir = self.plots_base_dir
    def setup_save_dirs(self):
        """
        Set up the directory structure:
        - self.base_dir: main data directory (self.save_dir)
        - self.plots_base_dir: main plots/logs directory
        - self.date_dir: date-stamped subdirectory for data
        - self.plots_date_dir: date-stamped subdirectory for plots/logs
        - self.plots_dir: 'plots' subdirectory inside plots_date_dir
        - self.md_path: path to the markdown log file (in plots_date_dir)
        """
        self.base_dir = self.save_dir
        self.plots_base_dir = self.plots_base_dir if hasattr(self, "plots_base_dir") else self.save_dir
        today = datetime.datetime.now().strftime("%Y-%m-%d")
        self.date_dir = os.path.join(self.base_dir, today)
        self.plots_date_dir = os.path.join(self.plots_base_dir, today)
        self.plots_dir = os.path.join(self.plots_date_dir, "plots")
        os.makedirs(self.date_dir, exist_ok=True)
        os.makedirs(self.plots_dir, exist_ok=True)
        self.md_path = os.path.join(self.plots_date_dir, f"{today}.md")
        if not os.path.exists(self.md_path):
            with open(self.md_path, "w") as f:
                f.write(f"# Simulation log for {today}\n\n")

    def log_file(self, filename):
        """
        Add a filename to the markdown log file.
        """
        if not hasattr(self, "md_path"):
            self.setup_save_dirs()
        with open(self.md_path, "a") as f:
            f.write(f"- {filename}\n")

    def save_plot_and_log(self, fig, plot_name, description=None):
        """
        Save a matplotlib figure to the plots directory and log it in the markdown file.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            The figure to save.
        plot_name : str
            The name of the plot file (without extension).
        description : str, optional
            Optional description to add to the markdown log.
        """
        if not hasattr(self, "plots_dir"):
            self.setup_save_dirs()
        plot_filename = self._get_unique_filename(plot_name, ".png", base_dir=self.plots_dir)
        print(f"Saving plot to {os.path.join(self.plots_dir, plot_filename)}")
        plot_path = os.path.join(self.plots_dir, plot_filename)
        fig.savefig(plot_path, bbox_inches="tight")
        rel_path = os.path.relpath(plot_path, self.plots_date_dir)
        with open(self.md_path, "a") as f:
            if description:
                f.write(f"\n### {plot_name}\n{description}\n")
            else:
                f.write(f"\n### {plot_name}\n")
            f.write(f"![{plot_name}]({rel_path})\n")
            f.write(f"- Saved as `{plot_filename}`\n")

    # Saving Data
    def _get_unique_filename(self, filename, ext, base_dir=None):
        """
        Generate a unique filename by appending a number if the file already exists.

        Parameters
        ----------
        filename : str
            Base filename (without extension).
        ext : str
            File extension (with dot), e.g., ".h5" or ".csv".

        Returns
        -------
        unique_filename : str
            Unique filename with extension.
        """
        if base_dir is None:
            base_dir = self.save_dir
        if not filename.endswith(ext):
            filename += ext
        base_filename, _ = os.path.splitext(filename)
        counter = 1
        unique_filename = filename
        while os.path.exists(os.path.join(base_dir, unique_filename)):
            unique_filename = f"{base_filename}_{counter}{ext}"
            counter += 1
        return unique_filename

    def _add_datetime_to_filename(self, filename, ext):
        """
        Add date and time to the filename before the extension.
        """
        now_str = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        base, _ = os.path.splitext(filename)
        return f"{base}_{now_str}{ext}"

    def save_results_to_hdf5(self, data_dict, filename=None):
        """
        Save arbitrary results to an HDF5 file.

        Parameters
        ----------
        data_dict : dict
            Dictionary where keys are dataset names and values are data (arrays, lists, etc.).
        filename : str, optional
            Name of the HDF5 file to save. Default is "results".
        """
        if filename is None:
            filename = self.sim_name
        if not hasattr(self, "save_dir"):
            raise AttributeError("The instance must have a 'save_dir' attribute specifying the save directory.")

        # Add date and time to filename
        filename = self._add_datetime_to_filename(filename, ".h5")
        filename = self._get_unique_filename(filename, ".h5", base_dir=self.save_dir)
        save_path = os.path.join(self.save_dir, filename)
        os.makedirs(self.save_dir, exist_ok=True)

        with h5py.File(save_path, "w") as f:
            for key, value in data_dict.items():
                arr = np.array([v.full() if hasattr(v, "full") else np.asarray(v) for v in value]) if isinstance(value, list) and value and (hasattr(value[0], "full") or hasattr(value[0], "__array__")) else np.asarray(value)
                f.create_dataset(key, data=arr)
        return filename

    def load_results_from_hdf5(self, filename="results", keys=None):
        """
        Load arbitrary results from an HDF5 file.

        Parameters
        ----------
        filename : str, optional
            Name of the HDF5 file to load. Default is "results.h5".
        keys : list of str, optional
            List of dataset names to load. If None, load all datasets.

        Returns
        -------
        results : dict
            Dictionary of loaded datasets.
        """
        if not filename.endswith(".h5"):
            filename += ".h5"
        file_path = os.path.join(self.save_dir, filename)
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"HDF5 file not found: {file_path}")

        results = {}
        with h5py.File(file_path, "r") as f:
            if keys is None:
                keys = list(f.keys())
            for key in keys:
                results[key] = np.array(f[key])
        return results

    def save_results_to_csv(self, data_dict, filename=None):
        """
        Save all results to a single CSV file using pandas.

        Parameters
        ----------
        data_dict : dict
            Dictionary where keys are dataset names and values are data (arrays, lists, etc.).
        filename : str, optional
            Name for the CSV file. Default is self.sim_name.
        """
        if filename is None:
            filename = self.sim_name
        if not hasattr(self, "save_dir"):
            raise AttributeError("The instance must have a 'save_dir' attribute specifying the save directory.")

        # Add date and time to filename
        filename = self._add_datetime_to_filename(filename, ".csv")
        filename = self._get_unique_filename(filename, ".csv", base_dir=self.save_dir)
        save_path = os.path.join(self.save_dir, filename)
        os.makedirs(self.save_dir, exist_ok=True)
        df = pd.DataFrame(data_dict)
        df.to_csv(save_path, index=False)
        print(f"Results saved to {save_path}")
        self.log_file(filename)
        return filename

    def append_to_csv(self, data_dict, filename=None):
        """
        Append new rows to an existing CSV file. Raises an error if the file does not exist.

        Parameters
        ----------
        data_dict : dict
            Dictionary where keys are column names and values are lists or scalars.
        filename : str, optional
            Name for the CSV file. Default is self.sim_name.
        """
        if filename is None:
            raise ValueError("filename cannot be None for append_to_csv.")
        if not filename.endswith(".csv"):
            filename += ".csv"
        save_path = os.path.join(self.save_dir, filename)
        os.makedirs(self.save_dir, exist_ok=True)

        if not os.path.exists(save_path):
            raise FileNotFoundError(f"CSV file not found: {save_path}")

        # Convert scalars to lists for consistency
        for k, v in data_dict.items():
            if not isinstance(v, (list, np.ndarray, pd.Series)):
                data_dict[k] = [v]

        df_new = pd.DataFrame(data_dict)
        df_existing = pd.read_csv(save_path)
        df_combined = pd.concat([df_existing, df_new], ignore_index=True)
        df_combined.to_csv(save_path, index=False)
        self.log_file(filename)

    def load_results_from_csv(self, filename):
        """
        Load all results from a single CSV file using pandas.

        Parameters
        ----------
        filename : str, optional
            Name for the CSV file. Default is self.sim_name.

        Returns
        -------
        results : dict
            Dictionary of loaded datasets.
        """
        if filename is None:
            raise ValueError("filename cannot be None for load_results_from_csv.")
        if not filename.endswith(".csv"):
            filename += ".csv"
        file_path = os.path.join(self.save_dir, filename)
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"CSV file not found: {file_path}")
        df = pd.read_csv(file_path)
        return df.to_dict(orient="list")

    def move_file(self, src_path, dest_dir):
        """
        Move a file from src_path to dest_dir.
        The file will be removed from the original directory.
        """
        if not os.path.isfile(src_path):
            raise FileNotFoundError(f"Source file not found: {src_path}")
        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir, exist_ok=True)
        dest_path = os.path.join(dest_dir, os.path.basename(src_path))
        shutil.move(src_path, dest_path)
        return dest_path
    
    def move_all_files(self, src_dir, dest_dir):
        """
        Move all files and folders from src_dir to dest_dir.
        """
        if not os.path.isdir(src_dir):
            raise NotADirectoryError(f"Source directory not found: {src_dir}")
        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir, exist_ok=True)
        for item in os.listdir(src_dir):
            src_path = os.path.join(src_dir, item)
            dest_path = os.path.join(dest_dir, item)
            if os.path.isfile(src_path):
                shutil.move(src_path, dest_path)
            elif os.path.isdir(src_path):
                shutil.move(src_path, dest_path)