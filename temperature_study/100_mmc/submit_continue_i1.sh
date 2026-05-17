#!/bin/bash

# Continue production for i1 from existing checkpoints.
TEMPS="295 300 305 310 315 320 325 330 335 340 345 350 355 360 365 370"

for T in $TEMPS; do
    DIR="${T}K/i1"
    DEFFNM="mmc_${T}K"

    if [ ! -d "$DIR" ]; then
        echo "[WARN] Missing directory: $DIR"
        continue
    fi

    if [ ! -f "$DIR/${DEFFNM}.tpr" ]; then
        echo "[WARN] Missing tpr: $DIR/${DEFFNM}.tpr"
        continue
    fi

    if [ ! -f "$DIR/${DEFFNM}.cpt" ]; then
        echo "[WARN] Missing checkpoint: $DIR/${DEFFNM}.cpt"
        continue
    fi

    echo "[RESUME] $DIR/$DEFFNM"
    (
        cd "$DIR" || exit 1
        sbatch --partition=standard --job-name="${DEFFNM}_cont" --time=96:00:00 --nodes=1 --ntasks-per-node=40 --cpus-per-task=1 \
               --mail-type=ALL --mail-user=user@institution.org \
               --wrap="module purge; module load GROMACS/2021.5-foss-2021b || module load GROMACS/2021.5 || module load GROMACS || module load gromacs; export GMX_MAXBACKUP=-1; export OMP_NUM_THREADS=1; srun --ntasks=40 gmx_mpi mdrun -deffnm ${DEFFNM} -cpi ${DEFFNM}.cpt -append -maxh 95"
    )
done
