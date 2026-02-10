#!/bin/bash

# Wrapper script to run run_anal_local.sh in the background using nohup.
# Usage: ./run_anal_local_bg.sh [iteration_number]
# Example: ./run_anal_local_bg.sh 1  -> Runs analysis in directory i1/ and logs to run_anal_local_i1.log

SCRIPT="run_anal_local.sh"

if [ ! -f "$SCRIPT" ]; then
    echo "Error: $SCRIPT not found in current directory."
    exit 1
fi

if [ -n "$1" ]; then
    ITER="$1"
    SUBDIR="i$ITER"
    LOGFILE="run_anal_local_${SUBDIR}.log"
    echo "Running analysis for iteration $ITER in subdirectory $SUBDIR"
    
    # Run with nohup, passing the subdirectory
    nohup bash "$SCRIPT" "$SUBDIR" > "$LOGFILE" 2>&1 &
else
    # Default behavior
    echo "No iteration number specified. Running in current directory."
    LOGFILE="run_anal_local.log"
    nohup bash "$SCRIPT" > "$LOGFILE" 2>&1 &
fi

echo "Starting $SCRIPT in the background..."
echo "Output will be redirected to $LOGFILE"

# Capture and display the Process ID (PID)
PID=$!
echo "Job started with PID: $PID"
echo "You can check the status with: ps -p $PID"
echo "You can follow the log with: tail -f $LOGFILE"
