#!/bin/bash

# Local workstation script for extending mmc_345K simulation
# Replicates the steps from extend.script

export GMX_MAXBACKUP=-1

# Define flags for local execution (consistent with run_local.sh)
MDRUN_FLAGS="-v -nt 32 -pin on --pinoffset 32"

echo "Starting local extension for mmc_345K..."
echo "MDRUN Flags: $MDRUN_FLAGS"

# 1) Extend the TPR file
echo "------------------------------------------------"
echo "Step 1: Extending TPR file to 500000 ps (500 ns)"
gmx convert-tpr -s mmc_345K.tpr -until 500000 -o mmc_345K.tpr

# 2) Continue the simulation
echo "------------------------------------------------"
echo "Step 2: Continuing simulation"
# -cpi specifies the checkpoint file to restart from
# -noappend creates new output files instead of appending
gmx mdrun $MDRUN_FLAGS -deffnm mmc_345K -cpi mmc_345K.cpt -noappend

echo "------------------------------------------------"
echo "Extension job complete."
