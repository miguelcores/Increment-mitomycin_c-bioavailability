#!/bin/bash

##----------------------- Start job description -----------------------

#SBATCH --partition=standard
#SBATCH --job-name=mmc_345K
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.cores@alumnos.upm.es

##----------------------- End job description ------------------------

# Function to check if running in Slurm
is_slurm() {
    if [ -n "$SLURM_JOB_ID" ]; then
        return 0
    else
        return 1
    fi
}

# Setup environment
if is_slurm; then
    echo "Detected Slurm environment. Loading modules..."
    module load apps/2021
    module load OpenMPI/4.1.4-GCC-12.2.0
    module load GROMACS/2021.5-foss-2021b
else
    echo "Running locally (no Slurm detected)."
fi

export GMX_MAXBACKUP=-1

# Define commands based on environment
if is_slurm; then
    # HPC / Slurm settings
    GMX_GROMPP="gmx grompp"
    GMX_MDRUN="mpirun gmx_mpi mdrun"
    MDRUN_FLAGS="-v"
else
    # Local workstation settings
    GMX_GROMPP="gmx grompp"
    GMX_MDRUN="gmx mdrun"
    # Flags from the commented out section in original script
    MDRUN_FLAGS="-v -nt 32 -pin on --pinoffset 32"
fi

echo "Using GMX_GROMPP: $GMX_GROMPP"
echo "Using GMX_MDRUN: $GMX_MDRUN"
echo "MDRUN Flags: $MDRUN_FLAGS"

# 1) Energy minimization
echo "------------------------------------------------"
echo "Step 1: Energy Minimization"
$GMX_GROMPP -f ../minim.mdp -c ../solvated.gro -p ../topol.top -o em.tpr -maxwarn 1
$GMX_MDRUN $MDRUN_FLAGS -deffnm em

# 2) NVT equilibration
echo "------------------------------------------------"
echo "Step 2: NVT Equilibration"
$GMX_GROMPP -f nvt.mdp -c em.gro -r em.gro -p ../topol.top -o nvt.tpr -maxwarn 1
$GMX_MDRUN $MDRUN_FLAGS -deffnm nvt

# 3) NPT equilibration
echo "------------------------------------------------"
echo "Step 3: NPT Equilibration"
$GMX_GROMPP -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p ../topol.top -o npt.tpr -maxwarn 1
$GMX_MDRUN $MDRUN_FLAGS -deffnm npt

# 4) Production at 345 K
echo "------------------------------------------------"
echo "Step 4: Production MD (345K)"
$GMX_GROMPP -f md.mdp -c npt.gro -t npt.cpt -p ../topol.top -o mmc_345K.tpr -maxwarn 1
$GMX_MDRUN $MDRUN_FLAGS -deffnm mmc_345K

echo "------------------------------------------------"
echo "Job complete."
