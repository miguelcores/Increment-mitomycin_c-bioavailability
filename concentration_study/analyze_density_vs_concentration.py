import argparse
import os
import re
from typing import List, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt


# Approximate atomic masses (amu or g/mol)
ATOMIC_MASS = {
    "H": 1.008,
    "C": 12.01,
    "N": 14.01,
    "O": 15.99,
    "S": 32.065,
    "P": 30.974,
    "Cl": 35.45,
}

# Water (SPC) mass: 2*H + O ≈ 18.015 amu
WATER_MASS = 18.015

# Avogadro's number
AVOGADRO = 6.02214076e23

# Conversion: 1 amu = 1.66054e-27 kg
AMU_TO_KG = 1.66054e-27


def parse_gro_box_volume(filepath: str) -> Optional[float]:
    """
    Extract box volume from .gro file (last line contains box vectors).
    Format: v1(x) v2(y) v3(z) [v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)]
    For cubic box: volume = v1x * v2y * v3z (converts nm^3 to m^3)
    """
    try:
        with open(filepath, "r") as f:
            lines = f.readlines()
            if len(lines) < 1:
                return None
            last_line = lines[-1].strip()
            parts = last_line.split()
            if len(parts) >= 3:
                v1x = float(parts[0])  # nm
                v2y = float(parts[1])  # nm
                v3z = float(parts[2])  # nm
                volume_nm3 = v1x * v2y * v3z
                volume_m3 = volume_nm3 * 1e-27  # Convert nm^3 to m^3
                return volume_m3
    except Exception:
        pass
    return None


def parse_topol_mass(filepath: str, n_mmc: int) -> Optional[float]:
    """
    Extract MMC mass from mmc.itp and water count from topol.top.
    Returns total mass in kg.
    """
    try:
        mmc_mass_amu = None
        n_water = None

        # Parse topol.top for molecule counts
        with open(filepath, "r") as f:
            in_molecules = False
            for line in f:
                line = line.strip()
                if line == "[ molecules ]":
                    in_molecules = True
                    continue
                if in_molecules:
                    if line.startswith("["):
                        break
                    if not line or line.startswith(";"):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        mol_name = parts[0]
                        mol_count = int(parts[1])
                        if mol_name.upper() == "SOL":
                            n_water = mol_count

        # Try to find MMC mass by parsing the included .itp file
        # Look for the mmc.itp include
        mmc_itp_path = None
        with open(filepath, "r") as f:
            for line in f:
                if "mmc.itp" in line.lower():
                    # Extract path from #include statement
                    match = re.search(r'#include\s+"([^"]+mmc\.itp)"', line)
                    if match:
                        rel_path = match.group(1)
                        mmc_itp_path = os.path.normpath(
                            os.path.join(os.path.dirname(filepath), rel_path)
                        )
                        break

        if mmc_itp_path and os.path.isfile(mmc_itp_path):
            mmc_mass_amu = sum_atom_masses_from_itp(mmc_itp_path)

        if mmc_mass_amu is None or n_water is None:
            return None

        # Convert to kg: each molecule mass in amu * AMU_TO_KG
        total_mass_kg = (n_mmc * mmc_mass_amu + n_water * WATER_MASS) * AMU_TO_KG
        return total_mass_kg

    except Exception:
        pass
    return None


def sum_atom_masses_from_itp(filepath: str) -> Optional[float]:
    """
    Sum atomic masses from [ atoms ] section in .itp file.
    """
    try:
        total_mass = 0.0
        in_atoms = False
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if line == "[ atoms ]":
                    in_atoms = True
                    continue
                if in_atoms:
                    if line.startswith("["):
                        break
                    if not line or line.startswith(";"):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        atom_type = parts[1]
                        # Try to extract element from atom type (e.g., "C" from "C1")
                        element = re.match(r"([A-Z][a-z]?)", atom_type)
                        if element:
                            elem = element.group(1)
                            if elem in ATOMIC_MASS:
                                total_mass += ATOMIC_MASS[elem]
        return total_mass if total_mass > 0 else None
    except Exception:
        pass
    return None


def parse_density_xvg(filepath: str) -> np.ndarray:
    data = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    time_ps = float(parts[0])
                    density = float(parts[1])
                    data.append((time_ps, density))
                except ValueError:
                    continue
    return np.array(data)


def find_density_files(root_dir: str) -> List[Tuple[int, str]]:
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
        density_path = os.path.join(full_path, "density.xvg")
        if os.path.isfile(density_path):
            found.append((concentration, density_path))

    return sorted(found, key=lambda x: x[0])


def compute_stats(density_data: np.ndarray, skip_fraction: float) -> Tuple[float, float, int]:
    if density_data.size == 0:
        return float("nan"), float("nan"), 0

    n = len(density_data)
    start_idx = int(n * skip_fraction)
    trimmed = density_data[start_idx:, 1] if n > 0 else np.array([])

    if trimmed.size == 0:
        return float("nan"), float("nan"), 0

    return float(np.mean(trimmed)), float(np.std(trimmed)), int(trimmed.size)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute system density vs MMC concentration from topol.top, solvated.gro, and density.xvg."
    )
    parser.add_argument(
        "--root",
        default=os.path.dirname(os.path.abspath(__file__)),
        help="Root directory containing *_mmc folders (default: script directory).",
    )
    parser.add_argument(
        "--skip-fraction",
        type=float,
        default=0.1,
        help="Fraction of initial frames to skip (default: 0.1).",
    )
    parser.add_argument(
        "--csv",
        default="density_vs_concentration.csv",
        help="Output CSV filename (default: density_vs_concentration.csv).",
    )
    parser.add_argument(
        "--plot",
        default="density_vs_concentration.png",
        help="Output plot filename (default: density_vs_concentration.png).",
    )

    args = parser.parse_args()
    root_dir = args.root

    density_files = find_density_files(root_dir)
    if not density_files:
        print("No density.xvg files found in *_mmc folders.")
        return

    results = []
    skipped = []

    for concentration, density_path in density_files:
        folder_path = os.path.dirname(density_path)
        topol_path = os.path.join(folder_path, "topol.top")
        solv_path = os.path.join(folder_path, "solvated.gro")

        # Get box volume
        volume_m3 = parse_gro_box_volume(solv_path)
        if volume_m3 is None:
            skipped.append((concentration, "No valid box volume"))
            continue

        # Get total mass
        total_mass_kg = parse_topol_mass(topol_path, concentration)
        if total_mass_kg is None:
            skipped.append((concentration, "Could not compute mass"))
            continue

        # Compute bulk density
        bulk_density_kg_m3 = total_mass_kg / volume_m3

        # Parse GROMACS-reported density
        data = parse_density_xvg(density_path)
        gmx_mean, gmx_std, gmx_n = compute_stats(data, args.skip_fraction)

        results.append(
            (
                concentration,
                bulk_density_kg_m3,
                gmx_mean,
                gmx_std,
                total_mass_kg,
                volume_m3,
                gmx_n,
                density_path,
            )
        )

    if not results:
        print("No valid results. Check for missing files or parsing errors.")
        return

    csv_path = os.path.join(root_dir, args.csv)
    with open(csv_path, "w", newline="") as f:
        f.write(
            "concentration_mmc,"
            "computed_bulk_density_kg_m3,"
            "gromacs_mean_density_kg_m3,"
            "gromacs_std_kg_m3,"
            "total_mass_kg,"
            "volume_m3,"
            "n_gromacs_points,"
            "source_file\n"
        )
        for row in results:
            f.write(
                f"{row[0]},"
                f"{row[1]:.6f},"
                f"{row[2]:.6f},"
                f"{row[3]:.6f},"
                f"{row[4]:.6f},"
                f"{row[5]:.9e},"
                f"{row[6]},"
                f"{os.path.relpath(row[7], root_dir)}\n"
            )

    concentrations = [r[0] for r in results]
    bulk_densities = [r[1] for r in results]
    gmx_means = [r[2] for r in results]
    gmx_stds = [r[3] for r in results]

    # Plot both computed bulk density and GROMACS-reported density
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.errorbar(
        concentrations,
        gmx_means,
        yerr=gmx_stds,
        fmt="o-",
        capsize=4,
        label="GROMACS (NPT output)",
        alpha=0.7,
    )
    ax.plot(
        concentrations,
        bulk_densities,
        "s--",
        label="Computed (total mass / volume)",
        alpha=0.7,
    )

    ax.axhline(y=1000, color="gray", linestyle=":", alpha=0.5, label="Pure water (~1000)")
    ax.set_xlabel("MMC concentration (molecule count)")
    ax.set_ylabel("Density (kg/m³)")
    ax.set_title("System density vs MMC concentration at 325 K")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    plot_path = os.path.join(root_dir, args.plot)
    plt.savefig(plot_path, dpi=150)
    plt.close()

    print(f"Wrote CSV: {csv_path}")
    print(f"Wrote plot: {plot_path}")

    if skipped:
        print("\nSkipped concentrations:")
        for conc, reason in skipped:
            print(f"  {conc}: {reason}")


if __name__ == "__main__":
    main()
