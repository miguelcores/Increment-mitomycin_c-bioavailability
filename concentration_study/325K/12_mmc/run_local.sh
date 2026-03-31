#!/bin/bash

# This script runs the simulation in the background using nohup.
# It effectively detaches the process so you can close the terminal.
# Logs are redirected to run.log

# Running with 32 threads (-nt 32) as requested.

nohup bash -c '
set -e

echo "Starting simulation for 12 MMC molecules..."

# 1) Energy minimization
echo "1. Energy Minimization"
gmx grompp -f ../mdp/minim.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 1
gmx mdrun -v -deffnm em -nt 32

# 2) NVT equilibration
echo "2. NVT Equilibration"
gmx grompp -f ../mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
gmx mdrun -deffnm nvt -nt 32

# 3) NPT equilibration
echo "3. NPT Equilibration"
gmx grompp -f ../mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
gmx mdrun -deffnm npt -nt 32

# 4) Production at 325 K
echo "4. Production Run"
gmx grompp -f ../mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -o mmc_12_325K.tpr -maxwarn 1
gmx mdrun -deffnm mmc_12_325K -nt 32

echo "Simulation 12_mmc complete."
' > run.log 2>&1 &

echo "Job 12_mmc submitted in background with 32 threads. Check run.log for output."
