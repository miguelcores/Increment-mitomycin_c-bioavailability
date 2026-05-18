#!/bin/bash
COUNTS="1 2 4 8 16 32 50 64 100"
for N in $(echo $COUNTS); do
    DIR="${N}_mmc/i0"
    SCRIPT="run_${N}_mmc.script"
    if [ ! -d "$DIR" ] || [ ! -f "$DIR/$SCRIPT" ]; then continue; fi
    echo "[SUBMIT] $DIR/$SCRIPT"
    (cd "$DIR" && sbatch "$SCRIPT")
done