#!/usr/bin/env bash
set -euo pipefail

# Usage: ./run_mmc_workflow.sh 32
NMMC="${1:-32}"
SEED="${2:-325}"

python3 prepare_mmc_system.py \
  --n-mmc "${NMMC}" \
  --seed "${SEED}" \
  --output-gro "system_${NMMC}mmc.gro" \
  --output-top "system_${NMMC}mmc.top" \
  --report "insertion_${NMMC}mmc.json"

python3 validate_mmc_system.py \
  --gro "system_${NMMC}mmc.gro" \
  --top "system_${NMMC}mmc.top" \
  --report "validation_${NMMC}mmc.json"

gmx grompp -f mdp/em_325K_1bar.mdp -c "system_${NMMC}mmc.gro" -p "system_${NMMC}mmc.top" -o em.tpr
gmx mdrun -deffnm em

gmx grompp -f mdp/nvt_325K_1bar.mdp -c em.gro -r em.gro -p "system_${NMMC}mmc.top" -o nvt.tpr
gmx mdrun -deffnm nvt

gmx grompp -f mdp/npt_325K_1bar.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p "system_${NMMC}mmc.top" -o npt.tpr
gmx mdrun -deffnm npt

gmx grompp -f mdp/md_325K_1bar.mdp -c npt.gro -t npt.cpt -p "system_${NMMC}mmc.top" -o md.tpr
gmx mdrun -deffnm md

# density check from trajectory/energies (interactive group selection piped)
printf "Density\n0\n" | gmx energy -f md.edr -o density_${NMMC}mmc.xvg
