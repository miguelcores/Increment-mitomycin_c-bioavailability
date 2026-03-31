#!/bin/bash

# Upload local temperature_study/100_mmc and submit all i1 jobs.

set -e

REMOTE_USER="${REMOTE_USER:-your_user}"
REMOTE_HOST="${REMOTE_HOST:-your.cluster.example.org}"
REMOTE_ROOT="${REMOTE_ROOT:-/path/to/remote/project/root}"
REMOTE_TS_DIR="$REMOTE_ROOT/temperature_study"
REMOTE_100_DIR="$REMOTE_TS_DIR/100_mmc"

LOCAL_ROOT="$(cd "$(dirname "$0")" && pwd)"
LOCAL_100_DIR="$LOCAL_ROOT/temperature_study/100_mmc"

if [ ! -d "$LOCAL_100_DIR" ]; then
    echo "Local folder not found: $LOCAL_100_DIR"
    exit 1
fi

CONTROL_PATH="/tmp/ssh-control-$USER-$$"
SSH_OPTS="-o ControlMaster=auto -o ControlPath=$CONTROL_PATH -o ControlPersist=10"

echo "Establishing SSH connection to $REMOTE_USER@$REMOTE_HOST ..."
ssh $SSH_OPTS -N -f "$REMOTE_USER@$REMOTE_HOST"

echo "Ensuring remote target directories exist..."
ssh $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST" "mkdir -p '$REMOTE_TS_DIR' '$REMOTE_100_DIR'"

echo "Uploading 100_mmc folder to $REMOTE_100_DIR ..."
ssh $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST" "rm -rf '$REMOTE_100_DIR' && mkdir -p '$REMOTE_100_DIR'"

# Stream-copy using tar over SSH to avoid scp canonicalization issues.
tar -C "$LOCAL_100_DIR" -cf - . | ssh $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST" "tar -xf - -C '$REMOTE_100_DIR'"

echo "Submitting i1 jobs remotely..."
ssh $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST" "REMOTE_ROOT='$REMOTE_ROOT' REMOTE_USER='$REMOTE_USER' bash -s" <<'REMOTE'
set -e
cd "${REMOTE_ROOT}/temperature_study/100_mmc"

if [ ! -f submit_all_i1.sh ]; then
    echo "submit_all_i1.sh not found in remote 100_mmc folder"
    exit 1
fi

echo "Running preflight grompp in 300K/i1 ..."
cd 300K/i1
module load apps/2021
module load OpenMPI/4.1.4-GCC-12.2.0
module load GROMACS/2021.5-foss-2021b
gmx grompp -f ../../minim.mdp -c ../../solvated.gro -p ../../topol.top -o preflight_em.tpr -maxwarn 1
rm -f preflight_em.tpr
cd ../../

sed -i 's/\r$//' submit_all_i1.sh
chmod +x submit_all_i1.sh
bash submit_all_i1.sh

echo "\nCurrent queue for $REMOTE_USER:"
squeue -u "$REMOTE_USER" || true

echo "\nMMC queue entries:"
squeue -u "$REMOTE_USER" | grep -E 'mmc_[0-9]+K' || true
REMOTE

echo "Closing SSH control connection..."
ssh -O exit -o ControlPath="$CONTROL_PATH" "$REMOTE_USER@$REMOTE_HOST" 2>/dev/null || true

echo "Done. Upload + i1 submission flow finished."
