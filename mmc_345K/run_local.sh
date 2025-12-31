#!/bin/bash

# Local workstation script for mmc_345K
# Replicates the steps from the Slurm script but for single-node/local execution

export GMX_MAXBACKUP=-1

# Define flags for local execution
# Adjust -nt (number of threads) and pinning as appropriate for your specific workstation
MDRUN_FLAGS="-v -nt 32 -pin on --pinoffset 32"

echo "Starting local execution for mmc_345K..."
echo "MDRUN Flags: $MDRUN_FLAGS"

# 1) Energy minimization
echo "------------------------------------------------"
echo "Step 1: Energy Minimization"
gmx grompp -f ../minim.mdp -c ../solvated.gro -p ../topol.top -o em.tpr -maxwarn 1
gmx mdrun $MDRUN_FLAGS -deffnm em

# 2) NVT equilibration
echo "------------------------------------------------"
echo "Step 2: NVT Equilibration"
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p ../topol.top -o nvt.tpr -maxwarn 1
gmx mdrun $MDRUN_FLAGS -deffnm nvt

# 3) NPT equilibration
echo "------------------------------------------------"
echo "Step 3: NPT Equilibration"
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p ../topol.top -o npt.tpr -maxwarn 1
gmx mdrun $MDRUN_FLAGS -deffnm npt

# 4) Production at 345 K
echo "------------------------------------------------"
echo "Step 4: Production MD (345K)"
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p ../topol.top -o mmc_345K.tpr -maxwarn 1
gmx mdrun $MDRUN_FLAGS -deffnm mmc_345K

echo "------------------------------------------------"
echo "Job complete."
