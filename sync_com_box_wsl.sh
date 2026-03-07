#!/bin/bash

# Configuration
REMOTE_USER="your_username"
REMOTE_HOST="your.cluster.hostname"
REMOTE_PATH="/path/to/remote/project"

# Example for alternative cluster:
# REMOTE_USER="alternative_user"
# REMOTE_HOST="alternative.cluster.hostname"
# REMOTE_PATH="/alternative/path"

LOCAL_ROOT="$(cd "$(dirname "$0")" && pwd)"

# Range configuration
MIN_CONC=12
MAX_CONC=150

# SSH connection sharing setup
CONTROL_PATH="/tmp/ssh-control-$USER-$$"
SSH_OPTS="-o ControlMaster=auto -o ControlPath=$CONTROL_PATH -o ControlPersist=10"

echo "Syncing com_all.xvg and box.xvg files from $REMOTE_HOST..."
echo "Range: ${MIN_CONC}_mmc to ${MAX_CONC}_mmc"
echo "You will be prompted for your password once."
echo ""

# Establish SSH connection (will prompt for password)
echo "Establishing SSH connection..."
ssh $SSH_OPTS -N -f "$REMOTE_USER@$REMOTE_HOST"
if [ $? -ne 0 ]; then
    echo "Failed to establish SSH connection. Exiting."
    exit 1
fi

echo "Connection established. Discovering remote concentration folders..."
mapfile -t folders < <(
    ssh $SSH_OPTS "$REMOTE_USER@$REMOTE_HOST" \
    "find '$REMOTE_PATH/concentration_study' -maxdepth 1 -type d -name '*_mmc' | sed 's#.*/##'" \
    | awk -v min="$MIN_CONC" -v max="$MAX_CONC" '
        /^[0-9]+_mmc$/ {
            split($0, a, "_")
            n = a[1] + 0
            if (n >= min && n <= max) print $0
        }
    ' | sort -V
)

if [ ${#folders[@]} -eq 0 ]; then
    echo "No remote folders found in range ${MIN_CONC}_mmc..${MAX_CONC}_mmc"
    ssh -O exit -o ControlPath="$CONTROL_PATH" "$REMOTE_USER@$REMOTE_HOST" 2>/dev/null
    exit 1
fi

echo "Found ${#folders[@]} folders to sync."
echo ""

com_pulled=0
box_pulled=0
missing_com=0
missing_box=0

for folder in "${folders[@]}"; do
    local_dir="$LOCAL_ROOT/concentration_study/$folder"
    remote_dir="$REMOTE_USER@$REMOTE_HOST:$REMOTE_PATH/concentration_study/$folder"

    mkdir -p "$local_dir"
    echo "Syncing $folder..."

    echo -n "  - com_all.xvg... "
    if scp $SSH_OPTS -q "$remote_dir/com_all.xvg" "$local_dir/" 2>/dev/null; then
        echo "✓"
        ((com_pulled++))
    else
        echo "✗ (missing or failed)"
        ((missing_com++))
    fi

    echo -n "  - box.xvg / boc.xvg... "
    if scp $SSH_OPTS -q "$remote_dir/box.xvg" "$local_dir/" 2>/dev/null; then
        echo "✓ (box.xvg)"
        ((box_pulled++))
    elif scp $SSH_OPTS -q "$remote_dir/boc.xvg" "$local_dir/" 2>/dev/null; then
        echo "✓ (boc.xvg)"
        ((box_pulled++))
    else
        echo "✗ (missing or failed)"
        ((missing_box++))
    fi
done

# Close SSH connection
ssh -O exit -o ControlPath="$CONTROL_PATH" "$REMOTE_USER@$REMOTE_HOST" 2>/dev/null

echo ""
echo "Sync complete."
echo "Folders processed: ${#folders[@]}"
echo "com_all.xvg: pulled=$com_pulled, missing=$missing_com"
echo "box/boc.xvg: pulled=$box_pulled, missing=$missing_box"
