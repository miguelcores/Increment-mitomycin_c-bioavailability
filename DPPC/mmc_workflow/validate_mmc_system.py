#!/usr/bin/env python3
"""Validate topology/coordinates after MMC insertion.

Checks:
1) molecule counts between GRO and [molecules]
2) total charge from molecule topology charges
3) estimated density from topology composition and box volume
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


@dataclass
class MoleculeDef:
    name: str
    charge: float
    mass: float


def parse_top_molecules(top_path: Path) -> Dict[str, int]:
    molecules: Dict[str, int] = {}
    in_mols = False
    for raw in top_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("["):
            in_mols = line.lower() == "[ molecules ]"
            continue
        if in_mols:
            parts = line.split()
            if len(parts) >= 2:
                molecules[parts[0]] = int(parts[1])
    if not molecules:
        raise RuntimeError(f"No [ molecules ] entries found in {top_path}")
    return molecules


def parse_molecule_from_itp(itp_path: Path, mol_name: str) -> MoleculeDef:
    in_atoms = False
    found_name = False
    total_q = 0.0
    total_m = 0.0

    for raw in itp_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("["):
            sec = line.lower()
            in_atoms = sec == "[ atoms ]" and found_name
            continue

        if not found_name:
            if re.match(rf"^{re.escape(mol_name)}\s+\d+", line):
                found_name = True
            continue

        if in_atoms:
            cols = line.split()
            if len(cols) < 8:
                continue
            total_q += float(cols[6])
            total_m += float(cols[7])

    if not found_name:
        raise RuntimeError(f"Moleculetype {mol_name} not found in {itp_path}")

    return MoleculeDef(mol_name, total_q, total_m)


def parse_gro_counts_and_box(gro_path: Path) -> Tuple[Dict[str, int], Tuple[float, float, float]]:
    lines = gro_path.read_text().splitlines()
    natoms = int(lines[1].strip())
    atom_lines = lines[2 : 2 + natoms]
    box_parts = lines[2 + natoms].split()
    if len(box_parts) < 3:
        raise RuntimeError("Only orthorhombic boxes are supported")
    box = (float(box_parts[0]), float(box_parts[1]), float(box_parts[2]))

    residues_seen = set()
    counts: Dict[str, int] = {}
    for l in atom_lines:
        resid = int(l[0:5])
        resname = l[5:10].strip()
        key = (resid, resname)
        if key in residues_seen:
            continue
        residues_seen.add(key)
        counts[resname] = counts.get(resname, 0) + 1

    return counts, box


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate DPPC+MMC inserted system")
    ap.add_argument("--gro", type=Path, default=Path("system_mmc.gro"))
    ap.add_argument("--top", type=Path, default=Path("system_mmc.top"))
    ap.add_argument("--dppc-itp", type=Path, default=Path("../dppc.itp"))
    ap.add_argument("--mmc-itp", type=Path, default=Path("mmc.itp"))
    ap.add_argument("--report", type=Path, default=Path("validation_report.json"))
    args = ap.parse_args()

    top_counts = parse_top_molecules(args.top)
    gro_counts, box = parse_gro_counts_and_box(args.gro)

    defs = {
        "DPPC": parse_molecule_from_itp(args.dppc_itp, "DPPC"),
        "MMC": parse_molecule_from_itp(args.mmc_itp, "MMC"),
        "SOL": MoleculeDef("SOL", 0.0, 18.0153),
    }

    comparison = {}
    for name in ("DPPC", "MMC", "SOL"):
        comparison[name] = {
            "top": top_counts.get(name, 0),
            "gro": gro_counts.get(name, 0),
            "match": top_counts.get(name, 0) == gro_counts.get(name, 0),
        }

    total_charge = 0.0
    total_mass_amu = 0.0
    for name, n in top_counts.items():
        if name not in defs:
            continue
        total_charge += defs[name].charge * n
        total_mass_amu += defs[name].mass * n

    volume_nm3 = box[0] * box[1] * box[2]
    volume_m3 = volume_nm3 * 1e-27
    mass_kg = total_mass_amu * 1.66053906660e-27
    density_kg_m3 = mass_kg / volume_m3

    result = {
        "counts": comparison,
        "top_charge_e": total_charge,
        "charge_ok": abs(total_charge) < 1e-6,
        "box_nm": list(box),
        "volume_nm3": volume_nm3,
        "estimated_density_kg_m3": density_kg_m3,
    }

    args.report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
