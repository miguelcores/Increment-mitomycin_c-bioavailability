#!/bin/bash

# List of temperatures to process (excluding 345)
TEMPS="300 305 310 315 320 325 330 335 340 350 355 360 365 370"

# Enable nullglob to handle case where directory might not exist
shopt -s nullglob

for T in $TEMPS; do
    # Temperature folders now use plain names: 300K, 305K, etc.
    MATCHES=(${T}K)
    
    # Check if we found a directory
    if [ ${#MATCHES[@]} -eq 0 ]; then
        echo "Warning: No directory found for temperature ${T}"
        continue
    fi
    
    # Take the first match (assuming only one valid folder per temp)
    DIR="${MATCHES[0]}"
    
    echo "Processing directory: $DIR"
    
    for SUB in i1 i2; do
        TARGET_DIR="$DIR/$SUB"
        
        if [ -d "$TARGET_DIR" ]; then
            # Script files keep the mmc_* naming convention.
            SCRIPT_NAME="mmc_${T}K_anal.script"
            
            echo "  Checking $TARGET_DIR for $SCRIPT_NAME..."
            
            if [ -f "$TARGET_DIR/$SCRIPT_NAME" ]; then
                echo "    Submitting $SCRIPT_NAME..."
                # Navigate to directory to execute sbatch so logs appear there
                (cd "$TARGET_DIR" && sbatch "$SCRIPT_NAME")
            else
                echo "    Error: Script $SCRIPT_NAME not found in $TARGET_DIR"
            fi
        else
            echo "  Skipping $SUB: Directory $TARGET_DIR does not exist"
        fi
    done
done
