#!/bin/bash

# Local workstation script for mmc_345K analysis
# Replicates the steps from mmc_345K_anal.script
# Usage: ./run_anal_local.sh [subdir]

# Check if a subdirectory argument is provided
if [ -n "$1" ]; then
    RUN_DIR="$1"
    echo "Analysis directory specified: $RUN_DIR"
    if [ -d "$RUN_DIR" ]; then
        cd "$RUN_DIR" || exit
        echo "Changed directory to $(pwd)"
    else
        echo "Error: Directory $RUN_DIR does not exist."
        exit 1
    fi
fi

export GMX_MAXBACKUP=-1

# Set number of threads for analysis tools (equivalent to cpus-per-task=32)
export OMP_NUM_THREADS=32

echo "Starting local analysis for mmc_345K..."
echo "GROMACS VERSION USED FOR ANALYSIS:"
gmx --version

# 1) Energy minimization: potential energy
echo "------------------------------------------------"
echo "Step 1: Energy Minimization Analysis"
# Example: if 'Potential' is entry 10 in gmx energy menu:
printf "10\n0\n" | gmx energy -f em.edr -o potential.xvg

# 2) NVT equilibration: temperature
echo "------------------------------------------------"
echo "Step 2: NVT Equilibration Analysis"
# Example: if 'Temperature' is entry 15:
printf "14\n0\n" | gmx energy -f nvt.edr -o temperature.xvg

# 3) NPT equilibration: pressure and density
echo "------------------------------------------------"
echo "Step 3: NPT Equilibration Analysis"
# Example: if 'Pressure' is 16 and 'Density' is 18:
printf "16\n0\n" | gmx energy -f npt.edr -o pressure.xvg
printf "22\n0\n" | gmx energy -f npt.edr -o density.xvg

# 4) MSD analysis of MMC (or X2IT group) from production
echo "------------------------------------------------"
echo "Step 4: MSD Analysis"
# First, create index file non-interactively.
# Example: if the group you want is number 2 in the default list, and you just want to keep it:
printf "2\nq\n" | gmx make_ndx -f mmc_345K.tpr -o index.ndx

# Now compute MSD using that group (assume it is group 2 in index.ndx)
# MSD: answer both prompts with group 2, then terminate selection with 0 or empty
printf "2\n" | gmx msd -f mmc_345K.xtc -s mmc_345K.tpr -n index.ndx -o msd_mmc_345K.xvg

echo "------------------------------------------------"
echo "Analysis complete."
