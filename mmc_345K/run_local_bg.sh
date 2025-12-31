#!/bin/bash

# Wrapper script to run run_local.sh in the background using nohup.
# This allows the job to continue running even if the terminal session is closed.

# Name of the script to run
SCRIPT="run_local.sh"
# Log file for output
LOGFILE="run_local.log"

if [ ! -f "$SCRIPT" ]; then
    echo "Error: $SCRIPT not found in current directory."
    exit 1
fi

echo "Starting $SCRIPT in the background..."
echo "Output will be redirected to $LOGFILE"

# Run with nohup
nohup bash "$SCRIPT" > "$LOGFILE" 2>&1 &

# Capture and display the Process ID (PID)
PID=$!
echo "Job started with PID: $PID"
echo "You can check the status with: ps -p $PID"
echo "You can follow the log with: tail -f $LOGFILE"
