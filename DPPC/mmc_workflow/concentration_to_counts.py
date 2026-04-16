#!/usr/bin/env python3
"""Convert target MMC concentration (mM) into molecule counts for this box."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

NA = 6.02214076e23


def parse_box_volume_nm3(gro: Path) -> float:
    lines = gro.read_text().splitlines()
    natoms = int(lines[1].strip())
    box = [float(x) for x in lines[2 + natoms].split()[:3]]
    return box[0] * box[1] * box[2]


def main() -> None:
    ap = argparse.ArgumentParser(description="Compute MMC molecule counts from target concentration")
    ap.add_argument("--gro", type=Path, default=Path("../dppc_325K_1bar.part0011.gro"))
    ap.add_argument("--dppc", type=int, default=256)
    ap.add_argument("--sol", type=int, default=7712)
    ap.add_argument("--waters-per-mmc", type=int, default=2)
    ap.add_argument("--mM", type=float, nargs="+", required=True, help="Target concentrations in mM")
    args = ap.parse_args()

    vol_nm3 = parse_box_volume_nm3(args.gro)
    vol_L = vol_nm3 * 1e-24

    rows = []
    for c_mM in args.mM:
        c_M = c_mM * 1e-3
        n_mmc = int(round(c_M * vol_L * NA))
        n_sol = args.sol - n_mmc * args.waters_per_mmc
        rows.append(
            {
                "target_mM": c_mM,
                "mmc_count": n_mmc,
                "dppc_count": args.dppc,
                "sol_count": n_sol,
                "valid": n_sol >= 0,
            }
        )

    print(json.dumps({"volume_nm3": vol_nm3, "counts": rows}, indent=2))


if __name__ == "__main__":
    main()
