#!/bin/bash

# Sync pre-extracted frame files for temperature_study/100_mmc/i0 and i1.

# Configuration
REMOTE_USER="HPC_USERNAME"
REMOTE_HOST="HPC_HOST"
REMOTE_ROOT="/home/w481/SCRATCH/miguelcores"
LOCAL_ROOT="$(cd "$(dirname "$0")" && pwd)"

# Must match the frame interval used when files were created remotely
FRAME_DT_PS=50000

# Temperatures and iterations to sync
temps=(295 300 305)
iterations=(i0)

# SSH connection sharing
CONTROL_PATH="/tmp/ssh-control-$USER-$$"
SSH_OPTS="-o ControlMaster=auto -o ControlPath=$CONTROL_PATH -o ControlPersist=10"

set -e

echo "Establishing SSH connection to $REMOTE_HOST..."
ssh $SSH_OPTS -N -f "$REMOTE_USER@$REMOTE_HOST"

echo "Syncing frame files to local machine..."
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

echo "Done. Synced files are under temperature_study/100_mmc/<T>K/<iX>/frames_mmc_<T>K_<iX>_${FRAME_DT_PS}ps.pdb"
