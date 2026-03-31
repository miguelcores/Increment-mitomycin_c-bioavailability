import argparse
import os
import re
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np


def read_xvg(filepath: str) -> np.ndarray:
    data = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    data.append([float(parts[0]), float(parts[1])])
                except ValueError:
                    continue
    return np.array(data)


def find_msd_files(root_dir: str, min_conc: int, max_conc: int) -> List[Tuple[int, str]]:
    found = []
    pattern = re.compile(r"^(\d+)_mmc$")

    for entry in os.listdir(root_dir):
        full_path = os.path.join(root_dir, entry)
        if not os.path.isdir(full_path):
            continue
        match = pattern.match(entry)
        if not match:
            continue
        concentration = int(match.group(1))
        if concentration < min_conc or concentration > max_conc:
            continue
        xvg_name = f"msd_mmc_{concentration}_325K.xvg"
        xvg_path = os.path.join(full_path, xvg_name)
        if os.path.isfile(xvg_path):
            found.append((concentration, xvg_path))

    return sorted(found, key=lambda x: x[0])


def compute_loglog_slope(
    time_ps: np.ndarray,
    msd_nm2: np.ndarray,
    fit_start_ns: float,
    fit_end_ns: float,
) -> Tuple[float, float, float, np.ndarray]:
    valid = (time_ps > 0) & (msd_nm2 > 0)
    time_ns = time_ps / 1000.0
    fit_mask = valid & (time_ns >= fit_start_ns) & (time_ns <= fit_end_ns)

    if np.count_nonzero(fit_mask) < 2:
        fit_mask = valid

    fit_time_ps = time_ps[fit_mask]
    fit_msd_nm2 = msd_nm2[fit_mask]

    log_time = np.log(fit_time_ps)
    log_msd = np.log(fit_msd_nm2)
    slope, intercept = np.polyfit(log_time, log_msd, 1)
    return slope, intercept, fit_time_ps.min(), fit_time_ps.max()


def plot_msd(concentration: int, filepath: str, fit_start_ns: float, fit_end_ns: float) -> None:
    data = read_xvg(filepath)
    if data.size == 0:
        print(f"No data found in {filepath}")
        return

    time_ps = data[:, 0]
    msd_nm2 = data[:, 1]

    output_dir = os.path.dirname(filepath)
    base_name = os.path.splitext(os.path.basename(filepath))[0]

    # Normal plot (time in ns)
    time_ns = time_ps / 1000.0
    plt.figure(figsize=(8, 6))
    plt.plot(time_ns, msd_nm2, color="blue", alpha=0.8, label=f"{concentration} MMC")
    plt.title(f"Mean Squared Displacement at {concentration} MMC")
    plt.xlabel("Time (ns)")
    plt.ylabel("MSD (nm^2)")
    plt.legend()
    plt.grid(True)

    normal_path = os.path.join(output_dir, f"{base_name}.png")
    plt.savefig(normal_path, dpi=150)
    plt.close()
    print(f"Saved plot to {normal_path}")

    # Log-log plot (time in ps)
    valid = (time_ps > 0) & (msd_nm2 > 0)
    if not np.any(valid):
        print(f"Skipping log-log plot for {filepath}: no positive values")
        return

    slope, intercept, fit_start_ps, fit_end_ps = compute_loglog_slope(
        time_ps, msd_nm2, fit_start_ns, fit_end_ns
    )
    fit_time = np.linspace(fit_start_ps, fit_end_ps, 200)
    fit_msd = np.exp(intercept) * (fit_time**slope)

    plt.figure(figsize=(8, 6))
    plt.loglog(time_ps[valid], msd_nm2[valid], ".", markersize=2, alpha=0.6)
    plt.loglog(
        fit_time,
        fit_msd,
        "r-",
        linewidth=2,
        label=f"alpha={slope:.3f} ({fit_start_ns:.0f}-{fit_end_ns:.0f} ns)",
    )
    plt.title(f"MSD Log-Log Plot at {concentration} MMC")
    plt.xlabel("Time (ps)")
    plt.ylabel("MSD (nm^2)")
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.5)

    log_path = os.path.join(output_dir, f"{base_name}_loglog.png")
    plt.savefig(log_path, dpi=150)
    plt.close()
    print(f"Saved log-log plot to {log_path}")


def plot_group_loglog(
    group_name: str,
    series: List[Tuple[int, str]],
    output_dir: str,
    fit_start_ns: float,
    fit_end_ns: float,
) -> None:
    if not series:
        return

    plt.figure(figsize=(12, 8))
    for concentration, filepath in series:
        data = read_xvg(filepath)
        if data.size == 0:
            continue
        time_ps = data[:, 0]
        msd_nm2 = data[:, 1]
        valid = (time_ps > 0) & (msd_nm2 > 0)
        if not np.any(valid):
            continue

        slope, intercept, fit_start_ps, fit_end_ps = compute_loglog_slope(
            time_ps, msd_nm2, fit_start_ns, fit_end_ns
        )
        fit_time = np.linspace(fit_start_ps, fit_end_ps, 200)
        fit_msd = np.exp(intercept) * (fit_time**slope)

        plt.loglog(time_ps[valid], msd_nm2[valid], alpha=0.7)
        plt.loglog(fit_time, fit_msd, linewidth=1.5, alpha=0.9)
        plt.plot([], [], label=f"mmc_{concentration} (alpha={slope:.3f})")

    plt.title(
        f"Log-Log Mean Square Displacement (alpha fit: {fit_start_ns:.0f}-{fit_end_ns:.0f} ns)"
    )
    plt.xlabel("Time (ps)")
    plt.ylabel("MSD (nm^2)")
    plt.grid(True, which="both", ls="-", alpha=0.4)
    plt.legend(fontsize=8)
    plt.tight_layout()

    output_path = os.path.join(output_dir, f"msd_loglog_{group_name}.png")
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Saved group plot to {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate normal and log-log MSD plots for concentration study."
    )
    parser.add_argument(
        "--root",
        default=os.path.dirname(os.path.abspath(__file__)),
        help="Root directory containing *_mmc folders (default: script directory).",
    )
    parser.add_argument("--min", dest="min_conc", type=int, default=2)
    parser.add_argument("--max", dest="max_conc", type=int, default=100)
    parser.add_argument("--fit-start-ns", type=float, default=50.0)
    parser.add_argument("--fit-end-ns", type=float, default=450.0)

    args = parser.parse_args()
    root_dir = args.root

    msd_files = find_msd_files(root_dir, args.min_conc, args.max_conc)
    if not msd_files:
        print("No MSD files found.")
        return

    for concentration, xvg_path in msd_files:
        plot_msd(concentration, xvg_path, args.fit_start_ns, args.fit_end_ns)

    group_low = [(c, p) for c, p in msd_files if 2 <= c <= 48]
    group_high = [(c, p) for c, p in msd_files if 52 <= c <= 100]
    plot_group_loglog(
        "2-48mmc",
        group_low,
        root_dir,
        args.fit_start_ns,
        args.fit_end_ns,
    )
    plot_group_loglog(
        "52-100mmc",
        group_high,
        root_dir,
        args.fit_start_ns,
        args.fit_end_ns,
    )


if __name__ == "__main__":
    main()
