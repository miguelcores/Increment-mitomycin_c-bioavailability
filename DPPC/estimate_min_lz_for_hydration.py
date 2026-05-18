#!/usr/bin/env python3
"""Estimate minimal box-height (Lz) increase to preserve hydration after MMC insertion.

This utility targets membrane systems where MMC is inserted by replacing water molecules
(e.g., `gmx insert-molecules -replace "resname SOL"`).

It uses:
- baseline DPPC/SOL counts from a base topology,
- baseline box vectors from a base GRO,
- optional empirical replacement-rate extraction from existing `*_mmc/i*/topol.top`.

Then it computes the minimum pre-insertion waters needed so that post-insertion waters
remain above `min_water_per_dppc * n_dppc`.
"""

from __future__ import annotations

import argparse
import glob
import math
import re
from pathlib import Path
from statistics import mean
from typing import Iterable


TOPOL_COUNT_RE = re.compile(r"^\s*([A-Za-z0-9_+-]+)\s+(\d+)\s*$")


def parse_topol_molecule_count(topol_path: Path, mol_name: str) -> int:
    count = None
    in_molecules = False

    for raw in topol_path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue

        if line.lower() == "[ molecules ]":
            in_molecules = True
            continue

        if not in_molecules:
            continue

        m = TOPOL_COUNT_RE.match(raw)
        if not m:
            continue

        name, value = m.group(1), int(m.group(2))
        if name == mol_name:
            count = value
            break

    if count is None:
        raise ValueError(f"Molecule '{mol_name}' not found in {topol_path}")

    return count


def parse_box_from_gro(gro_path: Path) -> tuple[float, float, float]:
    lines = gro_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid GRO file (too short): {gro_path}")

    box_tokens = lines[-1].split()
    if len(box_tokens) < 3:
        raise ValueError(f"Invalid GRO box line in {gro_path}: {lines[-1]!r}")

    lx, ly, lz = map(float, box_tokens[:3])
    return lx, ly, lz


def parse_headgroup_slab_geometry(
    gro_path: Path,
    residue_name: str,
    headgroup_atom_name: str,
) -> dict[str, float]:
    """Estimate the two water-slab thicknesses from headgroup planes.

    We split headgroup atoms into two leaflets along z (lower/upper halves).
    The lower and upper leaflet mean z-positions define planes that bound the
    hydrocarbon core. Water slabs are then from z=0 to lower-plane and from
    upper-plane to z=Lz.
    """

    lines = gro_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid GRO file (too short): {gro_path}")

    natoms = int(lines[1].strip())
    atoms = lines[2 : 2 + natoms]
    box_tokens = lines[2 + natoms].split()
    if len(box_tokens) < 3:
        raise ValueError(f"Invalid GRO box line in {gro_path}: {lines[2 + natoms]!r}")

    lx, ly, lz = map(float, box_tokens[:3])

    z_values: list[float] = []
    for atom in atoms:
        resname = atom[5:10].strip()
        atom_name = atom[10:15].strip()
        if resname == residue_name and atom_name == headgroup_atom_name:
            z_values.append(float(atom[36:44]))

    if len(z_values) < 4:
        raise ValueError(
            f"Not enough {headgroup_atom_name} atoms in residue {residue_name} to define leaflets"
        )

    z_values.sort()
    half = len(z_values) // 2
    lower = z_values[:half]
    upper = z_values[half:]

    lower_mean = mean(lower)
    upper_mean = mean(upper)

    water_bottom = max(0.0, lower_mean)
    water_top = max(0.0, lz - upper_mean)
    water_slab_height = water_bottom + water_top

    if water_slab_height <= 0.0:
        raise ValueError(
            "Computed non-positive total water-slab height from headgroup planes"
        )

    area = lx * ly

    return {
        "lx": lx,
        "ly": ly,
        "lz": lz,
        "area": area,
        "lower_plane_z": lower_mean,
        "upper_plane_z": upper_mean,
        "water_bottom": water_bottom,
        "water_top": water_top,
        "water_slab_height": water_slab_height,
    }


def extract_empirical_replacement_rates(
    base_sol: int,
    pattern: str,
) -> list[tuple[str, int, int, float]]:
    """Return tuples (path, n_mmc, replaced_waters, replaced_per_mmc)."""
    rows: list[tuple[str, int, int, float]] = []

    for top_path_str in sorted(glob.glob(pattern)):
        top_path = Path(top_path_str)

        mmc_match = re.search(r"(\d+)_mmc[/\\]i\d+[/\\]topol\.top$", str(top_path))
        if not mmc_match:
            continue

        n_mmc = int(mmc_match.group(1))
        if n_mmc <= 0:
            continue

        try:
            sol_count = parse_topol_molecule_count(top_path, "SOL")
            x2it_count = parse_topol_molecule_count(top_path, "X2IT")
        except ValueError:
            continue

        if x2it_count != n_mmc:
            continue

        replaced = base_sol - sol_count
        replaced_per_mmc = replaced / n_mmc
        rows.append((str(top_path), n_mmc, replaced, replaced_per_mmc))

    return rows


def compute_lz_targets(
    counts: Iterable[int],
    dppc_count: int,
    base_sol: int,
    base_lx: float,
    base_ly: float,
    base_lz: float,
    min_water_per_dppc: float,
    replace_per_mmc: float,
    water_number_density: float,
) -> list[dict[str, float]]:
    min_sol_post = math.ceil(dppc_count * min_water_per_dppc)
    area = base_lx * base_ly
    waters_per_nm_lz = area * water_number_density

    out: list[dict[str, float]] = []
    for n_mmc in counts:
        min_pre_sol = min_sol_post + replace_per_mmc * n_mmc
        delta_waters = min_pre_sol - base_sol
        delta_lz = max(0.0, delta_waters / waters_per_nm_lz)
        target_lz = base_lz + delta_lz

        out.append(
            {
                "n_mmc": float(n_mmc),
                "min_pre_sol": float(min_pre_sol),
                "delta_waters": float(delta_waters),
                "delta_lz": float(delta_lz),
                "target_lz": float(target_lz),
                "delta_each_side": float(delta_lz / 2.0),
            }
        )

    return out


def print_table(title: str, rows: list[dict[str, float]]) -> None:
    print(f"\n{title}")
    print("N_mmc  preSOL_target  dN_waters  dLz(nm)  Lz_target(nm)  +per_side(nm)")
    for row in rows:
        print(
            f"{int(row['n_mmc']):5d}"
            f"  {row['min_pre_sol']:13.1f}"
            f"  {row['delta_waters']:9.1f}"
            f"  {row['delta_lz']:7.4f}"
            f"  {row['target_lz']:13.4f}"
            f"  {row['delta_each_side']:12.4f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Estimate minimal box-height increase to maintain hydration threshold after MMC insertion."
    )
    parser.add_argument(
        "--base-topol",
        default="topol.top",
        help="Base topology path (default: topol.top)",
    )
    parser.add_argument(
        "--base-gro",
        default="dppc_325K_1bar.part0011.gro",
        help="Base GRO path (default: dppc_325K_1bar.part0011.gro)",
    )
    parser.add_argument(
        "--counts",
        nargs="+",
        type=int,
        default=[32, 50, 64, 100],
        help="MMC counts to evaluate (default: 32 50 64 100)",
    )
    parser.add_argument(
        "--min-water-per-dppc",
        type=float,
        default=30.0,
        help="Hydration target in waters per DPPC (default: 30)",
    )
    parser.add_argument(
        "--replace-per-mmc",
        type=float,
        default=None,
        help="Fixed waters replaced per MMC. If omitted, empirical mean is used when available.",
    )
    parser.add_argument(
        "--conservative-replace-per-mmc",
        type=float,
        default=None,
        help="Conservative waters replaced per MMC. If omitted, empirical max is used when available.",
    )
    parser.add_argument(
        "--water-number-density",
        type=float,
        default=33.3679,
        help="Water number density in molecules/nm^3 (default: 33.3679).",
    )
    parser.add_argument(
        "--density-model",
        choices=["slab-headgroup", "bulk"],
        default="slab-headgroup",
        help=(
            "How to estimate water density for dN/dLz conversion: "
            "'slab-headgroup' (default) uses two water slabs above/below bilayer, "
            "'bulk' uses --water-number-density directly."
        ),
    )
    parser.add_argument(
        "--headgroup-atom",
        default="P8",
        help="Headgroup atom name used to locate leaflet planes (default: P8).",
    )
    parser.add_argument(
        "--residue-name",
        default="DPPC",
        help="Lipid residue name in GRO file for leaflet detection (default: DPPC).",
    )
    parser.add_argument(
        "--empirical-pattern",
        default="*_mmc/i*/topol.top",
        help="Glob pattern for empirical MMC topologies (default: *_mmc/i*/topol.top)",
    )

    args = parser.parse_args()

    base_topol = Path(args.base_topol)
    base_gro = Path(args.base_gro)

    if not base_topol.exists():
        raise SystemExit(f"Base topology not found: {base_topol}")
    if not base_gro.exists():
        raise SystemExit(f"Base GRO not found: {base_gro}")

    dppc_count = parse_topol_molecule_count(base_topol, "DPPC")
    base_sol = parse_topol_molecule_count(base_topol, "SOL")
    lx, ly, lz = parse_box_from_gro(base_gro)

    slab_geom = parse_headgroup_slab_geometry(
        gro_path=base_gro,
        residue_name=args.residue_name,
        headgroup_atom_name=args.headgroup_atom,
    )

    slab_density = base_sol / (slab_geom["area"] * slab_geom["water_slab_height"])

    if args.density_model == "slab-headgroup":
        effective_water_density = slab_density
    else:
        effective_water_density = args.water_number_density

    empirical_rows = extract_empirical_replacement_rates(base_sol, args.empirical_pattern)

    empirical_mean = None
    empirical_max = None
    if empirical_rows:
        rates = [row[3] for row in empirical_rows]
        empirical_mean = mean(rates)
        empirical_max = max(rates)

    replace_per_mmc = args.replace_per_mmc
    if replace_per_mmc is None:
        replace_per_mmc = empirical_mean if empirical_mean is not None else 16.0

    conservative_replace_per_mmc = args.conservative_replace_per_mmc
    if conservative_replace_per_mmc is None:
        if empirical_max is not None:
            conservative_replace_per_mmc = empirical_max
        else:
            conservative_replace_per_mmc = replace_per_mmc + 1.5

    min_sol_post = math.ceil(dppc_count * args.min_water_per_dppc)
    area = lx * ly
    waters_per_nm_lz = area * effective_water_density

    print("Base system")
    print(f"  DPPC count: {dppc_count}")
    print(f"  Base SOL count: {base_sol}")
    print(f"  Base box: Lx={lx:.5f} Ly={ly:.5f} Lz={lz:.5f} nm")
    print(f"  Area Lx*Ly: {area:.4f} nm^2")
    print(
        "  Headgroup planes (from "
        f"{args.residue_name}:{args.headgroup_atom} means): "
        f"z_lower={slab_geom['lower_plane_z']:.4f}, z_upper={slab_geom['upper_plane_z']:.4f}"
    )
    print(
        "  Water slabs: "
        f"bottom={slab_geom['water_bottom']:.4f} nm, "
        f"top={slab_geom['water_top']:.4f} nm, "
        f"total={slab_geom['water_slab_height']:.4f} nm"
    )
    print(
        "  Slab-derived water density: "
        f"{slab_density:.4f} molecules/nm^3"
    )
    print(
        "  Density model used: "
        f"{args.density_model} (effective rho={effective_water_density:.4f})"
    )
    print(f"  Added waters per +1.0 nm Lz: {waters_per_nm_lz:.1f}")
    print(f"  Hydration target: >= {min_sol_post} SOL post-insertion")
    print(
        "  Slab-based relation: Nw = rho*A*(h_bottom + h_top), "
        "h_bottom' = h_bottom + dLz/2, h_top' = h_top + dLz/2, "
        "therefore dNw = rho*A*dLz"
    )

    if empirical_rows:
        print("\nEmpirical replacement data (from existing *_mmc/i*/topol.top)")
        print("path                               N_mmc  replaced  replaced_per_mmc")
        for path, n_mmc, replaced, replaced_per_mmc in empirical_rows:
            print(f"{path:34s} {n_mmc:5d} {replaced:9d} {replaced_per_mmc:17.3f}")
        print(
            f"Empirical replacement-rate mean={empirical_mean:.3f}, "
            f"max={empirical_max:.3f} waters/MMC"
        )

    rows_est = compute_lz_targets(
        counts=args.counts,
        dppc_count=dppc_count,
        base_sol=base_sol,
        base_lx=lx,
        base_ly=ly,
        base_lz=lz,
        min_water_per_dppc=args.min_water_per_dppc,
        replace_per_mmc=replace_per_mmc,
        water_number_density=effective_water_density,
    )

    rows_cons = compute_lz_targets(
        counts=args.counts,
        dppc_count=dppc_count,
        base_sol=base_sol,
        base_lx=lx,
        base_ly=ly,
        base_lz=lz,
        min_water_per_dppc=args.min_water_per_dppc,
        replace_per_mmc=conservative_replace_per_mmc,
        water_number_density=effective_water_density,
    )

    print_table(
        f"Estimated target Lz using replace_per_mmc={replace_per_mmc:.3f}",
        rows_est,
    )
    print_table(
        f"Conservative target Lz using replace_per_mmc={conservative_replace_per_mmc:.3f}",
        rows_cons,
    )

    print("\nSuggested command pattern per concentration (use conservative Lz if uncertain):")
    for row in rows_cons:
        n_mmc = int(row["n_mmc"])
        target_lz = row["target_lz"]
        print(
            "  "
            f"gmx editconf -f {base_gro} -o {n_mmc}_mmc/boxed.gro -c -box {lx:.5f} {ly:.5f} {target_lz:.4f}"
        )


if __name__ == "__main__":
    main()
