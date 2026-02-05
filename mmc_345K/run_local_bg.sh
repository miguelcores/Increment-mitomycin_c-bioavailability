#!/bin/bash

# Wrapper script to run run_local.sh in the background using nohup.
# Usage: ./run_local_bg.sh [subdir] [offset]
# Example: ./run_local_bg.sh i1 0   -> Runs in directory i1/ with offset 0
# Example: ./run_local_bg.sh i2 16  -> Runs in directory i2/ with offset 16

SCRIPT="run_local.sh"

if [ -n "$1" ]; then
    SUBDIR="$1"
    OFFSET="${2:-0}"
    LOGFILE="run_local_${SUBDIR}.log"
    
    echo "Running in subdirectory $SUBDIR with offset $OFFSET"
    
    # Pass subdirectory AND offset to run_local.sh
    nohup bash "$SCRIPT" "$SUBDIR" "$OFFSET" > "$LOGFILE" 2>&1 &
    
else
    # Default behavior (root directory of 345K)
    echo "No subdirectory specified. Running in current directory with offset 0."
    LOGFILE="run_local.log"
    nohup bash "$SCRIPT" "" "0" > "$LOGFILE" 2>&1 &
fi

if [ ! -f "$SCRIPT" ]; then
    echo "Error: $SCRIPT not found in current directory."
    exit 1
fi

echo "Starting $SCRIPT in the background..."
echo "Output will be redirected to $LOGFILE"

# Capture and display the Process ID (PID)
PID=$!
echo "Job started with PID: $PID"
echo "You can check the status with: ps -p $PID"
echo "You can follow the log with: tail -f $LOGFILE"
