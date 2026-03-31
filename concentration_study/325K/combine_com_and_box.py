#!/usr/bin/env python3
"""
Combine COM coordinates with box dimensions in the format expected for cluster analysis.

Output format:
  Line 1: Time (ps)
  Line 2: Box dimensions Lx Ly Lz (nm)
  Line 3-N+2: COM coordinates of each MMC molecule (x, y, z)
  [repeat for each frame]
"""

import sys
import subprocess
from pathlib import Path
import re

def extract_box_dimensions(tpr_file, xtc_file, output_file="box.xvg"):
    """Extract box dimensions from trajectory using gmx traj."""
    
    print(f"Extracting box dimensions...")
    
    cmd = [
        "gmx", "traj",
        "-f", str(xtc_file),
        "-s", str(tpr_file),
        "-ob", str(output_file)
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode == 0:
            print(f"✓ Box dimensions extracted: {output_file}")
            return True
        else:
            print(f"✗ Error extracting box: {result.stderr}")
            return False
    except Exception as e:
        print(f"✗ Error: {e}")
        return False

def parse_xvg_data(xvg_file):
    """Parse XVG file and return data lines (skip header)."""
    
    data = []
    with open(xvg_file, 'r') as f:
        for line in f:
            # Skip comments and header
            if line.startswith(('#', '@')):
                continue
            if line.strip():
                data.append(line.strip())
    return data

def combine_com_and_box(com_file, box_file, output_file, n_molecules=4):
    """
    Combine COM coordinates and box dimensions into row format.
    
    Parameters:
    -----------
    com_file : str
        XVG file with COM coordinates (Time X1 Y1 Z1 X2 Y2 Z2 ...)
    box_file : str
        XVG file with box dimensions (Time Lx Ly Lz ...)
    output_file : str
        Output file with combined data in row format
    n_molecules : int
        Number of MMC molecules
    """
    
    print(f"\nCombining COM coordinates and box dimensions...")
    print(f"  COM file: {com_file}")
    print(f"  Box file: {box_file}")
    print(f"  Output: {output_file}")
    
    # Parse both files
    com_data = parse_xvg_data(com_file)
    box_data = parse_xvg_data(box_file)
    
    if len(com_data) != len(box_data):
        print(f"Warning: Different number of frames in COM ({len(com_data)}) and box ({len(box_data)}) files")
        # Use minimum
        n_frames = min(len(com_data), len(box_data))
    else:
        n_frames = len(com_data)
    
    print(f"  Number of frames: {n_frames}")
    print(f"  Number of molecules: {n_molecules}")
    
    # Write combined output
    with open(output_file, 'w') as f:
        f.write("# Combined COM coordinates with box dimensions\n")
        f.write(f"# Format: Time | Box(Lx Ly Lz) | COM for each MMC (x y z)\n")
        f.write(f"# Number of molecules: {n_molecules}\n")
        f.write(f"# Each frame takes {n_molecules + 2} lines:\n")
        f.write(f"#   Line 1: Time (ps)\n")
        f.write(f"#   Line 2: Box dimensions Lx Ly Lz (nm)\n")
        f.write(f"#   Lines 3-{n_molecules+2}: COM coordinates of MMC i (x y z) (nm)\n\n")
        
        for i in range(n_frames):
            com_line = com_data[i].split()
            box_line = box_data[i].split()
            
            # Extract time, COM coordinates, and box dimensions
            time = float(com_line[0])
            
            # COM coordinates: columns 1-3, 4-6, 7-9, 10-12 for 4 molecules
            coms = []
            expected_cols = 1 + 3 * n_molecules
            if len(com_line) >= expected_cols:
                for mol_idx in range(n_molecules):
                    x = float(com_line[1 + 3*mol_idx])
                    y = float(com_line[1 + 3*mol_idx + 1])
                    z = float(com_line[1 + 3*mol_idx + 2])
                    coms.append((x, y, z))
            else:
                print(f"Warning: Frame {i} has {len(com_line)} columns, expected {expected_cols}")
                continue
            
            # Box dimensions: columns 1, 2, 3 (after time)
            if len(box_line) >= 4:
                box_x = float(box_line[1])
                box_y = float(box_line[2])
                box_z = float(box_line[3])
            else:
                print(f"Warning: Frame {i} box has {len(box_line)} columns, expected 4+")
                continue
            
            # Write in row format
            f.write(f"{time:.3f}\n")
            f.write(f"{box_x:.6f} {box_y:.6f} {box_z:.6f}\n")
            for x, y, z in coms:
                f.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
    
    print(f"\n✓ Combined file created: {output_file}")
    print(f"  Total frames written: {n_frames}")
    
    # Show sample
    print(f"\n--- Sample output (first frame) ---")
    with open(output_file, 'r') as f:
        line_count = 0
        for line in f:
            if not line.startswith('#'):
                print(line.rstrip())
                line_count += 1
                if line_count >= n_molecules + 2:
                    break

def process_system(system_dir, n_mmc):
    """Process one system: extract box and combine with COM."""
    
    system_path = Path(system_dir)
    tpr_file = system_path / f"mmc_{n_mmc}_325K.tpr"
    xtc_file = system_path / f"mmc_{n_mmc}_325K.xtc"
    com_file = system_path / "com_all.xvg"
    box_file = system_path / "box.xvg"
    output_file = system_path / "com_with_box.dat"
    
    print(f"\n{'='*70}")
    print(f"Processing: {system_dir} ({n_mmc} MMC molecules)")
    print(f"{'='*70}")
    
    # Check if COM file exists
    if not com_file.exists():
        print(f"✗ COM file not found: {com_file}")
        print(f"  Please generate it first using gmx traj with index file")
        return False
    
    # Check if box file exists
    if not box_file.exists():
        # Try to extract box dimensions if TPR/XTC available
        if tpr_file.exists() and xtc_file.exists():
            print(f"Attempting to extract box dimensions...")
            if not extract_box_dimensions(tpr_file, xtc_file, str(box_file)):
                return False
        else:
            print(f"✗ Box file not found: {box_file}")
            print(f"  Please run on HPC: bash extract_box_dimensions.sh {system_dir} {n_mmc}")
            print(f"  Then download box.xvg to {system_dir}/")
            return False
    else:
        print(f"✓ Using existing box file: {box_file}")
    
    # Combine COM and box
    combine_com_and_box(str(com_file), str(box_file), str(output_file), n_mmc)
    
    return True

def discover_systems(min_n=8, max_n=150, root=Path('.')):
    systems = []
    for child in root.iterdir():
        if not child.is_dir():
            continue
        match = re.fullmatch(r"(\d+)_mmc", child.name)
        if not match:
            continue
        n_mmc = int(match.group(1))
        if min_n <= n_mmc <= max_n:
            systems.append((child.name, n_mmc))
    systems.sort(key=lambda item: item[1])
    return systems

def main():
    """Process system(s)."""
    if len(sys.argv) > 1 and sys.argv[1] == "--range":
        if len(sys.argv) != 4:
            print("Usage: python combine_com_and_box.py --range <min_mmc> <max_mmc>")
            print("Example: python combine_com_and_box.py --range 8 150")
            sys.exit(1)

        min_n = int(sys.argv[2])
        max_n = int(sys.argv[3])
        systems = discover_systems(min_n=min_n, max_n=max_n)

        if not systems:
            print(f"No *_mmc folders found in range {min_n}-{max_n}")
            sys.exit(1)

        results = {}
        for system_name, n_mmc in systems:
            success = process_system(system_name, n_mmc)
            results[system_name] = "SUCCESS" if success else "FAILED"

        print(f"\n\n{'='*70}")
        print("SUMMARY")
        print(f"{'='*70}")
        for system, status in results.items():
            symbol = "✓" if status == "SUCCESS" else "✗"
            print(f"{symbol} {system}: {status}")

        failed = any(status == "FAILED" for status in results.values())
        sys.exit(1 if failed else 0)

    if len(sys.argv) > 1:
        # Single system
        if len(sys.argv) < 3:
            print("Usage: python combine_com_and_box.py <system_dir> <n_molecules>")
            print("Example: python combine_com_and_box.py 4_mmc 4")
            sys.exit(1)
        
        system_dir = sys.argv[1]
        n_mmc = int(sys.argv[2])
        success = process_system(system_dir, n_mmc)
        sys.exit(0 if success else 1)
    else:
        # Process default range
        systems = discover_systems(min_n=8, max_n=150)
        if not systems:
            print("No *_mmc folders found in default range 8-150")
            sys.exit(1)
        
        results = {}
        for system_name, n_mmc in systems:
            if Path(system_name).exists():
                success = process_system(system_name, n_mmc)
                results[system_name] = "SUCCESS" if success else "FAILED"
            else:
                print(f"\n⊘ {system_name} not found - skipping")
                results[system_name] = "SKIPPED"
        
        # Summary
        print(f"\n\n{'='*70}")
        print("SUMMARY")
        print(f"{'='*70}")
        for system, status in results.items():
            symbol = "✓" if status == "SUCCESS" else "✗" if status == "FAILED" else "⊘"
            print(f"{symbol} {system}: {status}")

if __name__ == "__main__":
    main()
