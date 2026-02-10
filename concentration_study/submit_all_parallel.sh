#!/bin/bash
echo "Submitting all jobs in parallel..."

echo "Submitting 2_mmc..."
cd 2_mmc
chmod +x run_local.sh
./run_local.sh
cd ..

echo "Submitting 4_mmc..."
cd 4_mmc
chmod +x run_local.sh
./run_local.sh
cd ..

echo "Submitting 8_mmc..."
cd 8_mmc
chmod +x run_local.sh
./run_local.sh
cd ..

echo "Submitting 12_mmc..."
cd 12_mmc
chmod +x run_local.sh
./run_local.sh
cd ..

echo "Submitting 16_mmc..."
cd 16_mmc
chmod +x run_local.sh
./run_local.sh
cd ..

