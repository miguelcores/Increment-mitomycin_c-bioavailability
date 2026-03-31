#!/bin/bash

# Submit all i0 production jobs for 100_mmc temperature study.
TEMPS="300 305 310 315 320 325 330 335 340 345 350"

for T in $TEMPS; do
    DIR="${T}K/i0"
    SCRIPT="mmc_${T}K.script"

    if [ ! -d "$DIR" ]; then
        echo "[WARN] Missing directory: $DIR"
        continue
    fi

    if [ ! -f "$DIR/$SCRIPT" ]; then
        echo "[WARN] Missing script: $DIR/$SCRIPT"
        continue
    fi

    echo "[SUBMIT] $DIR/$SCRIPT"
    (cd "$DIR" && sbatch "$SCRIPT")
done
