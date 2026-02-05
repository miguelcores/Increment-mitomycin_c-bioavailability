#!/bin/bash

# Local workstation script for mmc_345K
# Replicates the steps from the Slurm script but for single-node/local execution
# Usage: ./run_local.sh [subdir] [pin_offset]
# If [subdir] is provided (e.g., i1), it creates/uses that directory and runs the simulation there.
# [pin_offset] defaults to 0 if not provided.

export GMX_MAXBACKUP=-1

# Configuration
CORES_PER_JOB=48

# Paths to shared input files (relative to the run directory)
# Defaults for running in mmc_345K/
MINIM_MDP="../minim.mdp"
SOLVATED_GRO="../solvated.gro"
TOPOL_TOP="../topol.top"

# Check if a subdirectory argument is provided
if [ -n "$1" ]; then
    RUN_DIR="$1"
    echo "Run directory specified: $RUN_DIR"
    
    # Get offset from 2nd argument, default to 0
    PIN_OFFSET="${2:-0}"
    echo "Using Pin Offset: $PIN_OFFSET"

    if [ ! -d "$RUN_DIR" ]; then
        echo "Creating directory $RUN_DIR..."
        mkdir -p "$RUN_DIR"
    fi
    
    # Copy local config files to the subfolder if they don't exist there
    # We assume run_local.sh is being run from the parent dir containing these files
    for f in md.mdp nvt.mdp npt.mdp; do
        if [ ! -f "$RUN_DIR/$f" ]; then
            echo "Copying $f to $RUN_DIR/"
            cp "$f" "$RUN_DIR/"
        fi
    done
    
    # Change to the subfolder
    cd "$RUN_DIR" || exit
    
    # Update paths to point two levels up
    MINIM_MDP="../../minim.mdp"
    SOLVATED_GRO="../../solvated.gro"
    TOPOL_TOP="../../topol.top"
    
    echo "Changed directory to $(pwd)"
else
    echo "No subdirectory specified. Running in current directory."
    PIN_OFFSET="${2:-0}"
fi

# Define flags for local execution
MDRUN_FLAGS="-v -nt $CORES_PER_JOB -pin on --pinoffset $PIN_OFFSET"

echo "Starting local execution for mmc_345K..."
echo "MDRUN Flags: $MDRUN_FLAGS"
echo "Using Inputs: $MINIM_MDP, $SOLVATED_GRO, $TOPOL_TOP"

# 1) Energy minimization
echo "------------------------------------------------"
echo "Step 1: Energy Minimization"
gmx grompp -f "$MINIM_MDP" -c "$SOLVATED_GRO" -p "$TOPOL_TOP" -o em.tpr -maxwarn 1
gmx mdrun $MDRUN_FLAGS -deffnm em

# 2) NVT equilibration
echo "------------------------------------------------"
echo "Step 2: NVT Equilibration"
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p "$TOPOL_TOP" -o nvt.tpr -maxwarn 1
gmx mdrun $MDRUN_FLAGS -deffnm nvt

# 3) NPT equilibration
echo "------------------------------------------------"
echo "Step 3: NPT Equilibration"
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p "$TOPOL_TOP" -o npt.tpr -maxwarn 1
gmx mdrun $MDRUN_FLAGS -deffnm npt

# 4) Production at 345 K
echo "------------------------------------------------"
echo "Step 4: Production MD (345K)"
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p "$TOPOL_TOP" -o mmc_345K.tpr -maxwarn 1
gmx mdrun $MDRUN_FLAGS -deffnm mmc_345K

echo "------------------------------------------------"
echo "Job complete."
