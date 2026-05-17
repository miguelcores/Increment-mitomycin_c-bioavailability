#!/bin/bash

# Extract and sync trajectory frames for temperature_study/100_mmc/<iX>.

# Configuration
REMOTE_USER="HPC_USERNAME"
REMOTE_HOST="HPC_HOST"
REMOTE_ROOT="/home/r727/SCRATCH/miguelcores"
LOCAL_ROOT="$(cd "$(dirname "$0")" && pwd)"

# Frame extraction settings (ps)
# Two frames per 100 ns -> dt = 50 ns
FRAME_DT_PS=50000
OUTPUT_DIR="."

# Temperatures to process
temps=(315 320 325 330 335 340 345 350)
iterations=(i1)

# SSH connection sharing
CONTROL_PATH="/tmp/ssh-control-$USER-$$"
SSH_OPTS="-o ControlMaster=auto -o ControlPath=$CONTROL_PATH -o ControlPersist=10"

set -e

echo "Establishing SSH connection to $REMOTE_HOST..."
ssh $SSH_OPTS -N -f "$REMOTE_USER@$REMOTE_HOST"

echo "Running remote frame extraction for 100_mmc/<iX>..."
ssh $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST" bash -s <<'REMOTE'
set +e

module purge
module load GROMACS/2021.5-foss-2021b || module load GROMACS/2021.5 || module load GROMACS || module load gromacs

REMOTE_ROOT="/home/r727/SCRATCH/miguelcores"
FRAME_DT_PS=50000
OUTPUT_DIR="."
temps=(315 320 325 330 335 340 345 350)
iterations=(i1)

for t in "${temps[@]}"; do
  for it in "${iterations[@]}"; do
    folder="$REMOTE_ROOT/temperature_study/100_mmc/${t}K/${it}"
    xtc="$folder/mmc_${t}K.xtc"
    tpr="$folder/mmc_${t}K.tpr"

    if [ ! -f "$xtc" ] || [ ! -f "$tpr" ]; then
      echo "[SKIP] ${t}K/${it}: missing .xtc or .tpr"
      continue
    fi

    mkdir -p "$folder/$OUTPUT_DIR"

    out="$folder/$OUTPUT_DIR/frames_mmc_${t}K_${it}_${FRAME_DT_PS}ps.pdb"
    # Extract a single multi-model PDB with snapshots every FRAME_DT_PS
    printf "0\n0\n" | gmx trjconv -s "$tpr" -f "$xtc" -o "$out" -dt "$FRAME_DT_PS" -pbc mol -center -ur compact >/dev/null 2>&1
    if [ $? -ne 0 ]; then
      echo "[SKIP] ${t}K/${it}: trjconv failed"
      continue
    fi

    echo "[OK] ${t}K/${it}: wrote frames_mmc_${t}K_${it}_${FRAME_DT_PS}ps.pdb"
  done
done
REMOTE

echo "Syncing extracted frames to local machine..."
for t in "${temps[@]}"; do
  for it in "${iterations[@]}"; do
    local_dir="$LOCAL_ROOT/temperature_study/100_mmc/${t}K/${it}"
    remote_file="$REMOTE_ROOT/temperature_study/100_mmc/${t}K/${it}/frames_mmc_${t}K_${it}_${FRAME_DT_PS}ps.pdb"

    mkdir -p "$local_dir"

    if ssh $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST" "[ -f '$remote_file' ]"; then
      scp -q $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST:$remote_file" "$local_dir/" || true
      echo "[SYNC] ${t}K/${it}"
    else
      echo "[SKIP] ${t}K/${it}: no frames file"
    fi
  done
done

# Close SSH connection
ssh -O exit -o ControlPath="$CONTROL_PATH" "$REMOTE_USER@$REMOTE_HOST" 2>/dev/null

echo "Done. Frames are under temperature_study/100_mmc/<T>K/<iX>/frames_mmc_<T>K_<iX>_${FRAME_DT_PS}ps.pdb"
