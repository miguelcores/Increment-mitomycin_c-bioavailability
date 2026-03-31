#!/bin/bash
set -e

for d in 8_mmc 16_mmc 32_mmc; do
  echo "Submitting ${d} at 300K..."
  cd "$d"
  sbatch run_${d%%_*}_mmc.script
  cd ..
done