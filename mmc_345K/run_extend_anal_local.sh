#!/bin/bash

# Local workstation script for mmc_345K extension analysis
# Replicates the steps from extend_anal.script

export GMX_MAXBACKUP=-1

# Set number of threads for analysis tools
export OMP_NUM_THREADS=32

echo "Starting local extension analysis for mmc_345K..."
echo "GROMACS VERSION USED FOR ANALYSIS:"
gmx --version

echo "------------------------------------------------"
echo "Step 1: Concatenate trajectories"
# creates mmc_345K_0-500ns.xtc containing both the initial run and the 500 ns extension
gmx trjcat -f mmc_345K.xtc mmc_345K.part0002.xtc -o mmc_345K_0-500ns.xtc

echo "------------------------------------------------"
echo "Step 2: Compute MSD"
# Now compute MSD using that group (assume it is group 2 in index.ndx)
# MSD: answer both prompts with group 2, then terminate selection with 0 or empty
printf "2\n" | gmx msd -f mmc_345K_0-500ns.xtc -s mmc_345K.tpr -n index.ndx -o msd_mmc_345K_0-500ns.xvg

echo "Analysis complete."
