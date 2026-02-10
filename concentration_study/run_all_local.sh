#!/bin/bash
set -e

echo "Entering 2_mmc..."
cd 2_mmc
chmod +x run_local.sh
./run_local.sh
cd ..

echo "Entering 4_mmc..."
cd 4_mmc
chmod +x run_local.sh
./run_local.sh
cd ..

echo "Entering 8_mmc..."
cd 8_mmc
chmod +x run_local.sh
./run_local.sh
cd ..

echo "Entering 12_mmc..."
cd 12_mmc
chmod +x run_local.sh
./run_local.sh
cd ..

echo "Entering 16_mmc..."
cd 16_mmc
chmod +x run_local.sh
./run_local.sh
cd ..

