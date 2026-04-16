#!/usr/bin/env python3
"""Prepare DPPC+MMC systems by replacing solvent water with MMC under overlap constraints.

This workflow keeps box dimensions and lipid count unchanged.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


@dataclass
class Atom:
    resid: int
    resname: str
    atomname: str
    atomnr: int
    x: float
    y: float
    z: float
    vx: float | None = None
    vy: float | None = None
    vz: float | None = None


def parse_gro(path: Path) -> tuple[str, List[Atom], tuple[float, float, float], bool]:
    lines = path.read_text().splitlines()
    title = lines[0].rstrip("\n")
    natoms = int(lines[1].strip())
    atom_lines = lines[2 : 2 + natoms]
    box_parts = lines[2 + natoms].split()
    if len(box_parts) < 3:
        raise ValueError("Only orthorhombic boxes are supported (3 box vectors expected).")
    box = (float(box_parts[0]), float(box_parts[1]), float(box_parts[2]))

    atoms: List[Atom] = []
    has_velocities = False
    for line in atom_lines:
        resid = int(line[0:5])
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomnr = int(line[15:20])
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        vx = vy = vz = None
        if len(line) >= 68:
            has_velocities = True
            vx = float(line[44:52])
            vy = float(line[52:60])
            vz = float(line[60:68])
        atoms.append(Atom(resid, resname, atomname, atomnr, x, y, z, vx, vy, vz))

    return title, atoms, box, has_velocities


def write_gro(path: Path, title: str, atoms: Sequence[Atom], box: tuple[float, float, float], write_velocities: bool) -> None:
    out = [title, f"{len(atoms):5d}"]
    for i, a in enumerate(atoms, start=1):
        resid = a.resid % 100000
        atomnr = i % 100000
        if write_velocities:
            vx = 0.0 if a.vx is None else a.vx
            vy = 0.0 if a.vy is None else a.vy
            vz = 0.0 if a.vz is None else a.vz
            out.append(
                f"{resid:5d}{a.resname:<5}{a.atomname:>5}{atomnr:5d}"
                f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{vx:8.4f}{vy:8.4f}{vz:8.4f}"
            )
        else:
            out.append(
                f"{resid:5d}{a.resname:<5}{a.atomname:>5}{atomnr:5d}"
                f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}"
            )
    out.append(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}")
    path.write_text("\n".join(out) + "\n")


def group_residues(atoms: Sequence[Atom]) -> Dict[tuple[int, str], List[int]]:
    groups: Dict[tuple[int, str], List[int]] = {}
    for i, a in enumerate(atoms):
        key = (a.resid, a.resname)
        groups.setdefault(key, []).append(i)
    return groups


def pbc_delta(dx: float, box: float) -> float:
    return dx - box * round(dx / box)


def min_image_dist2(a: tuple[float, float, float], b: tuple[float, float, float], box: tuple[float, float, float]) -> float:
    dx = pbc_delta(a[0] - b[0], box[0])
    dy = pbc_delta(a[1] - b[1], box[1])
    dz = pbc_delta(a[2] - b[2], box[2])
    return dx * dx + dy * dy + dz * dz


def random_rotation_matrix(rng: random.Random) -> List[List[float]]:
    u1 = rng.random()
    u2 = rng.random()
    u3 = rng.random()
    q1 = math.sqrt(1.0 - u1) * math.sin(2.0 * math.pi * u2)
    q2 = math.sqrt(1.0 - u1) * math.cos(2.0 * math.pi * u2)
    q3 = math.sqrt(u1) * math.sin(2.0 * math.pi * u3)
    q4 = math.sqrt(u1) * math.cos(2.0 * math.pi * u3)

    return [
        [1 - 2 * (q3 * q3 + q4 * q4), 2 * (q2 * q3 - q1 * q4), 2 * (q2 * q4 + q1 * q3)],
        [2 * (q2 * q3 + q1 * q4), 1 - 2 * (q2 * q2 + q4 * q4), 2 * (q3 * q4 - q1 * q2)],
        [2 * (q2 * q4 - q1 * q3), 2 * (q3 * q4 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3)],
    ]


def apply_rotation(v: tuple[float, float, float], r: Sequence[Sequence[float]]) -> tuple[float, float, float]:
    return (
        r[0][0] * v[0] + r[0][1] * v[1] + r[0][2] * v[2],
        r[1][0] * v[0] + r[1][1] * v[1] + r[1][2] * v[2],
        r[2][0] * v[0] + r[2][1] * v[1] + r[2][2] * v[2],
    )


def center_of_geometry(coords: Sequence[tuple[float, float, float]]) -> tuple[float, float, float]:
    n = len(coords)
    return (
        sum(c[0] for c in coords) / n,
        sum(c[1] for c in coords) / n,
        sum(c[2] for c in coords) / n,
    )


def wrap_to_box(v: tuple[float, float, float], box: tuple[float, float, float]) -> tuple[float, float, float]:
    return (v[0] % box[0], v[1] % box[1], v[2] % box[2])


def make_cell_key(c: tuple[float, float, float], box: tuple[float, float, float], ncell: tuple[int, int, int]) -> tuple[int, int, int]:
    ix = int((c[0] / box[0]) * ncell[0]) % ncell[0]
    iy = int((c[1] / box[1]) * ncell[1]) % ncell[1]
    iz = int((c[2] / box[2]) * ncell[2]) % ncell[2]
    return (ix, iy, iz)


def build_spatial_grid(
    coords_by_idx: Sequence[tuple[float, float, float]],
    indices: Sequence[int],
    box: tuple[float, float, float],
    cutoff: float,
) -> tuple[Dict[tuple[int, int, int], List[int]], tuple[int, int, int]]:
    ncell = (
        max(1, int(box[0] / cutoff)),
        max(1, int(box[1] / cutoff)),
        max(1, int(box[2] / cutoff)),
    )
    grid: Dict[tuple[int, int, int], List[int]] = {}
    for idx in indices:
        key = make_cell_key(coords_by_idx[idx], box, ncell)
        grid.setdefault(key, []).append(idx)
    return grid, ncell


def has_overlap_with_grid(
    trial_coords: Sequence[tuple[float, float, float]],
    coords_by_idx: Sequence[tuple[float, float, float]],
    grid: Dict[tuple[int, int, int], List[int]],
    ncell: tuple[int, int, int],
    box: tuple[float, float, float],
    min_dist2: float,
    ignore_indices: set[int],
) -> bool:
    for tc in trial_coords:
        cx, cy, cz = make_cell_key(tc, box, ncell)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    nk = (
                        (cx + dx) % ncell[0],
                        (cy + dy) % ncell[1],
                        (cz + dz) % ncell[2],
                    )
                    for idx in grid.get(nk, []):
                        if idx in ignore_indices:
                            continue
                        if min_image_dist2(tc, coords_by_idx[idx], box) < min_dist2:
                            return True
    return False


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Replace water molecules with MMC in DPPC bilayer system.")
    p.add_argument("--input-gro", type=Path, default=Path("../dppc_325K_1bar.part0011.gro"))
    p.add_argument("--mmc-gro", type=Path, default=Path("../../mitomycin_c/mmc.gro"))
    p.add_argument("--top-template", type=Path, default=Path("topol_mmc.top"))
    p.add_argument("--output-gro", type=Path, default=Path("system_mmc.gro"))
    p.add_argument("--output-top", type=Path, default=Path("system_mmc.top"))
    p.add_argument("--report", type=Path, default=Path("insertion_report.json"))
    p.add_argument("--n-mmc", type=int, required=True, help="Number of MMC molecules to insert")
    p.add_argument("--waters-per-mmc", type=int, default=2, help="Number of water molecules removed per MMC")
    p.add_argument("--min-dist", type=float, default=0.22, help="Minimum allowed distance (nm) between inserted MMC atoms and existing atoms")
    p.add_argument("--max-anchor-tries", type=int, default=2000)
    p.add_argument("--max-rot-tries", type=int, default=200)
    p.add_argument("--seed", type=int, default=325)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)

    title, atoms, box, has_vel = parse_gro(args.input_gro)
    _, mmc_atoms_raw, _, _ = parse_gro(args.mmc_gro)

    residue_groups = group_residues(atoms)
    dpcc_res = [k for k in residue_groups if k[1] == "DPPC"]
    water_res = [k for k in residue_groups if k[1] == "SOL"]

    if len(dpcc_res) != 256:
        raise RuntimeError(f"Expected 256 DPPC molecules, found {len(dpcc_res)}")

    waters_needed = args.n_mmc * args.waters_per_mmc
    if waters_needed > len(water_res):
        raise RuntimeError(
            f"Not enough water molecules: need {waters_needed}, available {len(water_res)}"
        )

    # Build a water lookup with oxygen position as insertion anchor.
    water_ow: Dict[tuple[int, str], tuple[float, float, float]] = {}
    for key in water_res:
        idxs = residue_groups[key]
        ow = [atoms[i] for i in idxs if atoms[i].atomname == "OW"]
        if not ow:
            continue
        water_ow[key] = (ow[0].x, ow[0].y, ow[0].z)

    if len(water_ow) != len(water_res):
        raise RuntimeError("Some SOL residues do not contain OW atoms.")

    # MMC template centered at COG.
    mmc_template_names = [a.atomname for a in mmc_atoms_raw]
    mmc_template_coords = [(a.x, a.y, a.z) for a in mmc_atoms_raw]
    cog = center_of_geometry(mmc_template_coords)
    mmc_centered = [(c[0] - cog[0], c[1] - cog[1], c[2] - cog[2]) for c in mmc_template_coords]

    kept_atom_indices = set(range(len(atoms)))
    removed_waters: List[tuple[int, str]] = []
    inserted_mmc_coords: List[List[tuple[float, float, float]]] = []

    available_waters = list(water_res)

    def nearest_waters(anchor_key: tuple[int, str]) -> List[tuple[int, str]]:
        anchor = water_ow[anchor_key]
        ranked = sorted(
            available_waters,
            key=lambda w: min_image_dist2(anchor, water_ow[w], box),
        )
        return ranked[: args.waters_per_mmc]

    existing_positions = [(a.x, a.y, a.z) for a in atoms]
    min_dist2 = args.min_dist * args.min_dist

    for ins_id in range(args.n_mmc):
        placed = False
        grid, ncell = build_spatial_grid(existing_positions, sorted(kept_atom_indices), box, args.min_dist)
        for _ in range(args.max_anchor_tries):
            if not available_waters:
                break
            anchor_key = rng.choice(available_waters)
            selected_waters = nearest_waters(anchor_key)

            ignore_indices: set[int] = set()
            for wkey in selected_waters:
                for idx in residue_groups[wkey]:
                    ignore_indices.add(idx)

            anchor = water_ow[anchor_key]

            for _rot in range(args.max_rot_tries):
                rot = random_rotation_matrix(rng)
                trial = []
                for c in mmc_centered:
                    rc = apply_rotation(c, rot)
                    pc = wrap_to_box((rc[0] + anchor[0], rc[1] + anchor[1], rc[2] + anchor[2]), box)
                    trial.append(pc)

                bad = has_overlap_with_grid(
                    trial,
                    existing_positions,
                    grid,
                    ncell,
                    box,
                    min_dist2,
                    ignore_indices,
                )

                if not bad:
                    for mmc_coords in inserted_mmc_coords:
                        for tc in trial:
                            for ec in mmc_coords:
                                if min_image_dist2(tc, ec, box) < min_dist2:
                                    bad = True
                                    break
                            if bad:
                                break
                        if bad:
                            break

                if not bad:
                    for wkey in selected_waters:
                        if wkey in available_waters:
                            available_waters.remove(wkey)
                        if wkey not in removed_waters:
                            removed_waters.append(wkey)
                        for idx in residue_groups[wkey]:
                            kept_atom_indices.discard(idx)
                    inserted_mmc_coords.append(trial)
                    placed = True
                    break

            if placed:
                break

        if not placed:
            raise RuntimeError(
                f"Could not place MMC #{ins_id + 1} without overlap. Try increasing --waters-per-mmc or lowering --min-dist."
            )

    kept_atoms = [atoms[i] for i in sorted(kept_atom_indices)]

    # Rebuild residue IDs with new MMC residues appended.
    max_resid = max(a.resid for a in kept_atoms)
    out_atoms: List[Atom] = []
    for a in kept_atoms:
        out_atoms.append(Atom(a.resid, a.resname, a.atomname, 0, a.x, a.y, a.z, a.vx, a.vy, a.vz))

    for i, mmc_coords in enumerate(inserted_mmc_coords, start=1):
        resid = max_resid + i
        for aname, c in zip(mmc_template_names, mmc_coords):
            out_atoms.append(Atom(resid, "MMC", aname, 0, c[0], c[1], c[2], 0.0 if has_vel else None, 0.0 if has_vel else None, 0.0 if has_vel else None))

    write_gro(args.output_gro, f"{title} + {args.n_mmc} MMC", out_atoms, box, has_vel)

    n_sol_final = len([k for k in residue_groups if k[1] == "SOL"]) - len(removed_waters)
    top_text = args.top_template.read_text()
    top_text = top_text.replace("__MMC_COUNT__", str(args.n_mmc))
    top_text = top_text.replace("__SOL_COUNT__", str(n_sol_final))
    args.output_top.write_text(top_text)

    report = {
        "input_gro": str(args.input_gro),
        "output_gro": str(args.output_gro),
        "output_top": str(args.output_top),
        "n_mmc_inserted": args.n_mmc,
        "waters_removed": len(removed_waters),
        "waters_per_mmc": args.waters_per_mmc,
        "final_sol": n_sol_final,
        "lipids_kept": len(dpcc_res),
        "box_nm": list(box),
        "min_dist_nm": args.min_dist,
        "seed": args.seed,
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n")

    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
