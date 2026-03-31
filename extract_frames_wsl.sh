#!/bin/bash

# Configuration
REMOTE_USER="${REMOTE_USER:-your_user}"
REMOTE_HOST="${REMOTE_HOST:-your.cluster.example.org}"
REMOTE_ROOT="${REMOTE_ROOT:-/path/to/remote/project/root}"
LOCAL_ROOT="$(cd "$(dirname "$0")" && pwd)"

# Frame extraction settings (ps)
# Two frames per 100 ns -> dt = 50 ns
FRAME_DT_PS=50000
OUTPUT_DIR="."

# Concentrations to process
concentrations=(2 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84 88 92 96 100)

# SSH connection sharing
CONTROL_PATH="/tmp/ssh-control-$USER-$$"
SSH_OPTS="-o ControlMaster=auto -o ControlPath=$CONTROL_PATH -o ControlPersist=10"

set -e

echo "Establishing SSH connection to $REMOTE_HOST..."
ssh $SSH_OPTS -N -f "$REMOTE_USER@$REMOTE_HOST"

echo "Running remote frame extraction..."
ssh $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST" bash -s <<'REMOTE'
set +e

module purge
module load GROMACS/2021.5-foss-2021b || module load GROMACS/2021.5 || module load GROMACS || module load gromacs

REMOTE_ROOT="/path/to/remote/project/root"
FRAME_DT_PS=50000
OUTPUT_DIR="."

concentrations=(2 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84 88 92 96 100)

for conc in "${concentrations[@]}"; do
    folder="$REMOTE_ROOT/concentration_study/${conc}_mmc"
    xtc="$folder/mmc_${conc}_325K.xtc"
    tpr="$folder/mmc_${conc}_325K.tpr"

    if [ ! -f "$xtc" ] || [ ! -f "$tpr" ]; then
        echo "[SKIP] ${conc}_mmc: missing .xtc or .tpr"
        continue
    fi

    mkdir -p "$folder/$OUTPUT_DIR"

    out="$folder/$OUTPUT_DIR/frames_mmc_${conc}_325K_${FRAME_DT_PS}ps.pdb"
    # Extract a single multi-model PDB with snapshots every FRAME_DT_PS
    # Provide both selection prompts (e.g., "System" twice)
    printf "0\n0\n" | gmx trjconv -s "$tpr" -f "$xtc" -o "$out" -dt "$FRAME_DT_PS" -pbc mol -center -ur compact >/dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "[SKIP] ${conc}_mmc: trjconv failed"
        continue
    fi

    echo "[OK] ${conc}_mmc: wrote frames_mmc_${conc}_325K_${FRAME_DT_PS}ps.pdb"
done
REMOTE

echo "Syncing frames to local machine..."
for conc in "${concentrations[@]}"; do
    folder="${conc}_mmc"
    local_dir="$LOCAL_ROOT/concentration_study/$folder"
    remote_file="$REMOTE_ROOT/concentration_study/$folder/frames_mmc_${conc}_325K_${FRAME_DT_PS}ps.pdb"

    mkdir -p "$local_dir"

    if ssh $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST" "[ -f '$remote_file' ]"; then
        scp -q $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST:$remote_file" "$local_dir/" || true
    else
        echo "[SKIP] ${folder}: no frames file"
    fi
done

# Close SSH connection
ssh -O exit -o ControlPath="$CONTROL_PATH" "$REMOTE_USER@$REMOTE_HOST" 2>/dev/null

echo "Done. Frames are under concentration_study/<N>_mmc/frames_mmc_<N>_325K_${FRAME_DT_PS}ps.pdb"