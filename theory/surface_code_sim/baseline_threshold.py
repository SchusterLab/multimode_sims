"""
Baseline 2D Rotated Surface Code Threshold Sweep
Reproduces Figure 11 (top-left) from Duckering et al. 2020
"Virtualized Logical Qubits: A 2.5D Architecture for Error-Corrected Quantum Computing"

Target threshold: p_th ≈ 0.009
"""

import stim
import sinter
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def make_surface_code_circuit(d: int, p: float) -> stim.Circuit:
    """
    Rotated surface code circuit using Stim's built-in generator.
    Error model: symmetric depolarizing (matches paper baseline).
    """
    return stim.Circuit.generated(
        "surface_code:unrotated_memory_z",
        rounds=d,
        distance=d,
        after_clifford_depolarization=p,
        after_reset_flip_probability=p,
        before_measure_flip_probability=p,
        before_round_data_depolarization=p,
    )


def run_threshold_sweep(
    distances: list[int] = [3, 5, 7],
    num_p_points: int = 10,
    p_min: float = 1e-3,
    p_max: float = 0.05,
    max_shots: int = 10_000,
    max_errors: int = 100,
    num_workers: int = 4,
    save_dir: str = "results",
) -> list[sinter.TaskStats]:
    """
    Run threshold sweep using sinter for parallelized sampling + pymatching decoding.
    """
    error_rates = np.logspace(np.log10(p_min), np.log10(p_max), num_p_points)

    tasks = [
        sinter.Task(
            circuit=make_surface_code_circuit(d, p),
            json_metadata={"d": d, "p": p},
        )
        for d in distances
        for p in error_rates
    ]

    print(f"Running {len(tasks)} tasks ({len(distances)} distances × {num_p_points} error rates)")
    print(f"max_shots={max_shots:,}, max_errors={max_errors}, num_workers={num_workers}")

    results = sinter.collect(
        tasks=tasks,
        num_workers=num_workers,
        max_shots=max_shots,
        max_errors=max_errors,
        decoders=["pymatching"],
        print_progress=True,
        save_resume_filepath=f"{save_dir}/checkpoint.csv",
    )

    return results


def plot_threshold(results: list[sinter.TaskStats], save_path: str = "threshold_baseline.png"):
    fig, ax = plt.subplots(figsize=(7, 5))

    sinter.plot_error_rate(
        ax=ax,
        stats=results,
        x_func=lambda stat: stat.json_metadata["p"],
        group_func=lambda stat: f"d={stat.json_metadata['d']}",
        failure_units_per_shot_func=lambda stat: 1,
    )

    ax.axvline(x=0.009, color="k", linestyle="--", alpha=0.4, label="p_th ≈ 0.009 (paper)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Physical Error Rate p", fontsize=12)
    ax.set_ylabel("Logical Error Rate", fontsize=12)
    ax.set_title("2D Rotated Surface Code — Baseline Threshold\n(cf. Duckering et al. 2020, Fig. 11)", fontsize=11)
    ax.legend(title="Distance")
    ax.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    print(f"Saved: {save_path}")
    plt.show()


def quick_check(distances: list[int] = [3, 5, 7], p: float = 0.005, shots: int = 10_000):
    """
    Quick sanity check using manual decode loop (no sinter).
    Below threshold (p=0.005 < 0.009): larger d should give lower logical error rate.
    """
    import pymatching

    print(f"\nQuick check at p={p} (below threshold, expect larger d → lower error rate):")
    for d in distances:
        circuit = make_surface_code_circuit(d, p)
        dem = circuit.detector_error_model(decompose_errors=True)
        matcher = pymatching.Matching.from_detector_error_model(dem)

        sampler = circuit.compile_detector_sampler()
        detection_events, obs_flips = sampler.sample(shots, separate_observables=True)
        predictions = matcher.decode_batch(detection_events)
        rate = np.mean(predictions != obs_flips)
        print(f"  d={d}, p={p} → logical error rate = {rate:.5f}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Surface code threshold sweep")
    parser.add_argument("--quick", action="store_true", help="Run quick sanity check only")
    parser.add_argument("--max-shots", type=int, default=10_000)
    parser.add_argument("--max-errors", type=int, default=500)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    save_dir = Path("results")
    save_dir.mkdir(exist_ok=True)

    if args.quick:
        quick_check()
    else:
        quick_check()
        results = run_threshold_sweep(
            max_shots=args.max_shots,
            max_errors=args.max_errors,
            num_workers=args.workers,
            save_dir=str(save_dir),
        )
        plot_threshold(results, save_path=str(save_dir / "threshold_baseline.png"))
