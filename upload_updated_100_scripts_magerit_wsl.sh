#!/bin/bash

# Upload only updated simulation scripts (and key control files) to Magerit,
# preserving folder structure under temperature_study/100_mmc.

set -e

REMOTE_USER="${REMOTE_USER:-your_user}"
REMOTE_HOST="${REMOTE_HOST:-your.cluster.example.org}"
REMOTE_ROOT="${REMOTE_ROOT:-}"

LOCAL_ROOT="$(cd "$(dirname "$0")" && pwd)"
LOCAL_100_DIR="$LOCAL_ROOT/temperature_study/100_mmc"

if [ ! -d "$LOCAL_100_DIR" ]; then
    echo "Local folder not found: $LOCAL_100_DIR"
    exit 1
fi

CONTROL_PATH="/tmp/ssh-control-$USER-$$"
SSH_OPTS="-o ControlMaster=auto -o ControlPath=$CONTROL_PATH -o ControlPersist=10"

cleanup() {
    ssh -O exit -o ControlPath="$CONTROL_PATH" "$REMOTE_USER@$REMOTE_HOST" 2>/dev/null || true
    [ -n "${TMP_LIST:-}" ] && [ -f "$TMP_LIST" ] && rm -f "$TMP_LIST"
}
trap cleanup EXIT

echo "Establishing SSH connection to $REMOTE_USER@$REMOTE_HOST ..."
ssh $SSH_OPTS -N -f "$REMOTE_USER@$REMOTE_HOST"

if [ -z "$REMOTE_ROOT" ]; then
    echo "Detecting remote root..."
    REMOTE_ROOT=$(ssh $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST" "
if [ -d \"\$HOME/SCRATCH/project_root\" ]; then
    printf '%s' \"\$HOME/SCRATCH/project_root\"
elif [ -d \"\$HOME/project_root\" ]; then
    printf '%s' \"\$HOME/project_root\"
else
    printf '%s' \"\$HOME/remote_project_root\"
fi
")
fi

REMOTE_100_DIR="$REMOTE_ROOT/temperature_study/100_mmc"

echo "Remote root: $REMOTE_ROOT"
echo "Remote target: $REMOTE_100_DIR"

# Build relative file list to upload (ONLY .script files in temp/iteration dirs).
TMP_LIST=$(mktemp)
(
    cd "$LOCAL_100_DIR"
    find . -type f \
    \( -path "./*K/i*/mmc_*K.script" \
        -o -path "./*K/i*/mmc_*K_anal.script" \) \
        | sed 's#^./##' | sort
) > "$TMP_LIST"

FILE_COUNT=$(wc -l < "$TMP_LIST" | tr -d ' ')
if [ "$FILE_COUNT" -eq 0 ]; then
    echo "No files matched upload rules. Nothing to upload."
    exit 0
fi

echo "Preparing remote directory..."
ssh $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST" "mkdir -p '$REMOTE_100_DIR'"

echo "Uploading $FILE_COUNT files (structure-preserving) ..."
# Stream tar with only selected files.
tar -C "$LOCAL_100_DIR" -cf - -T "$TMP_LIST" | ssh $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST" "tar -xf - -C '$REMOTE_100_DIR'"

echo "Normalizing line endings on remote scripts..."
ssh $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST" "find '$REMOTE_100_DIR' -type f \( -name '*.sh' -o -name '*.script' \) -exec sed -i 's/\r$//' {} +"

echo "Upload complete. Updated files are in $REMOTE_100_DIR"
