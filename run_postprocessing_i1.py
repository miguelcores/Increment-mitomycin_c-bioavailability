import os
import glob
import re
import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Import analysis functions from existing scripts
# We assume they are in the same directory
import analyze_msd_log
import generate_msd_plots

def get_i1_file(root_dir, temp):
    """
    Finds the msd xvg file for a given temperature in the i1 directory.
    Handles 'K' vs 'k' case sensitivity in folder names.
    Returns path if found, else None.
    """
    # Possible folder names: mmc_300K, mmc_300k, mmc_300, etc.
    # We specifically want the one ending in K/k for this project structure based on listing
    candidates = [
        f"mmc_{temp}",
        f"mmc_{temp.lower()}", # handle 315K vs 315k
    ]
    
    # Try to find the directory from candidates (and partial matching if needed)
    # The file list showed 'mmc_315k', 'mmc_300K'.
    # We can glob for the directory.
    
    # Construct glob pattern for directory: mmc_300[Kk]*
    temp_num = temp.replace('K', '')
    pattern = os.path.join(root_dir, f"mmc_{temp_num}[Kk]*")
    dirs = glob.glob(pattern)
    
    if not dirs:
        return None
    
    # Prefer exact match if possible, or take the first one
    # We need subdirectory "i1"
    for d in dirs:
        i1_dir = os.path.join(d, "i1")
        if os.path.isdir(i1_dir):
            # Look for xvg file
            # Naming convention: msd_mmc_300K.xvg
            # But script used to have _0-500ns. Check for any msd_mmc_*.xvg
            xvg_pattern = os.path.join(i1_dir, "msd_mmc_*.xvg")
            xvgs = glob.glob(xvg_pattern)
            if xvgs:
                # If multiple, prefer the one matching temp exactly if possible
                for xvg in xvgs:
                    if f"_{temp}." in xvg or f"_{temp}_" in xvg:
                         return xvg
                return xvgs[0] # Return first if no strict match
    
    return None

def plot_arrhenius_i1(csv_files):
    """
    Adapted from plot_arrhenius.py to work with a list of specific CSV files.
    """
    temps = []
    g_ds = []
    g_errs = []
    c_ds = []
    c_errs = []
    
    print("\nExtracting data for Arrhenius plot...")
    
    # Sort files by temperature extracted from filename or content
    # We'll rely on reading the CSV content which has "Temperature" field
    
    data_list = []
    
    for filepath in csv_files:
        try:
            with open(filepath, 'r') as f:
                reader = csv.reader(f)
                file_data = {}
                for row in reader:
                    if len(row) == 2:
                        file_data[row[0]] = row[1]
                
                temp_str = file_data.get('Temperature', 'N/A')
                if temp_str == 'N/A':
                    continue
                    
                # Clean temp string (remove K)
                try:
                    t_val = float(re.sub(r'[^\d\.]', '', temp_str))
                except:
                    continue
                
                g_d_val, g_d_err = parse_d_string(file_data.get('GROMACS D (1e-5 cm^2/s)', 'N/A'))
                c_d_val, c_d_err = parse_d_string(file_data.get('Calculated D (1e-5 cm^2/s)', 'N/A'))
                
                if c_d_val is not None:
                    data_list.append({
                        'T': t_val,
                        'G_D': g_d_val, 'G_Err': g_d_err,
                        'C_D': c_d_val, 'C_Err': c_d_err
                    })
        except Exception as e:
            print(f"Error reading CSV {filepath}: {e}")

    # Sort by T
    data_list.sort(key=lambda x: x['T'])
    
    # Prepare data for plotting
    calc_data = [] # (temp, d, err)
    grom_data = [] # (temp, d, err)

    for d in data_list:
        print(f"{d['T']}K: Calc D = {d['C_D']} +/- {d['C_Err']}")
        # Always have Calc D based on previous check
        calc_data.append((d['T'], d['C_D'], d['C_Err']))
        
        if d['G_D'] is not None:
             grom_data.append((d['T'], d['G_D'], d['G_Err']))

    if not calc_data:
        print("No valid data for Arrhenius plot.")
        return

    # Create arrays for Calculated D
    c_temps = np.array([x[0] for x in calc_data])
    c_inv_temps = 1000.0 / c_temps
    c_ds = np.array([x[1] for x in calc_data])
    c_errs = np.array([x[2] for x in calc_data])
    ln_c_ds = np.log(c_ds)
    ln_c_errs = c_errs / c_ds

    # Create arrays for GROMACS D
    if grom_data:
        g_temps = np.array([x[0] for x in grom_data])
        g_inv_temps = 1000.0 / g_temps
        g_ds = np.array([x[1] for x in grom_data])
        g_errs = np.array([x[2] for x in grom_data])
        ln_g_ds = np.log(g_ds)
        ln_g_errs = g_errs / g_ds
    else:
        g_inv_temps = None
    
    # Plotting side-by-side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: GROMACS D
    if grom_data:
        ax1.errorbar(g_inv_temps, ln_g_ds, yerr=ln_g_errs, fmt='o', label='GROMACS D', capsize=5, color='blue')
        
        # Fit
        slope, intercept, r_value, p_value, std_err = stats.linregress(g_inv_temps, ln_g_ds)
        fit_line = slope * g_inv_temps + intercept
        ax1.plot(g_inv_temps, fit_line, 'r--', label=f'Fit: R²={r_value**2:.4f}\nSlope={slope:.4f}')
        
        ax1.set_title('Arrhenius Plot (i1 Runs): GROMACS D')
        ax1.set_xlabel('1000 / T (K⁻¹)')
        ax1.set_ylabel('ln(D)')
        ax1.legend()
        ax1.grid(True)
    else:
        ax1.text(0.5, 0.5, 'No GROMACS D Data', ha='center')

    # Plot 2: Calculated D
    ax2.errorbar(c_inv_temps, ln_c_ds, yerr=ln_c_errs, fmt='o', label='Calculated D', capsize=5, color='green')
    
    result = stats.linregress(c_inv_temps, ln_c_ds)
    slope = result.slope
    intercept = result.intercept
    r_value = result.rvalue
    
    fit_line = slope * c_inv_temps + intercept
    ax2.plot(c_inv_temps, fit_line, 'r--', label=f'Fit: R²={r_value**2:.4f}\nSlope={slope:.4f}')
    
    ax2.set_title('Arrhenius Plot (i1 Runs): Calculated D')
    ax2.set_xlabel('1000 / T (K⁻¹)')
    ax2.set_ylabel('ln(D)')
    ax2.legend()
    ax2.grid(True)
    
    plt.tight_layout()
    output_path = 'arrhenius_plots_i1.png'
    plt.savefig(output_path)
    plt.close()
    print(f"\nSaved Arrhenius plot to {output_path}")

def parse_d_string(d_str):
    # Format: "0.1234 +/- 0.01" or "N/A"
    if d_str == 'N/A' or d_str is None:
        return None, None
    try:
        parts = d_str.split('+/-')
        val = float(parts[0].strip())
        err = float(parts[1].strip()) if len(parts) > 1 else 0.0
        return val, err
    except:
        return None, None

def main():
    root_dir = os.getcwd()
    # Temperatures from 295 to 370, step 5 (excluding 345? user said "from 295k to 370k" but noted 345 was special/skipped in scripts? previous prompt: "but mmc_345". I'll skip 345 if requested, but analysis should probably include it if exists. The user said "but mmc_345" in context of creating submission script, here "Do the postprocessing... from 295k to 370k". I will include 345 if found.)
    # Actually, previous prompt "make script that goes and executes sbatch mmc_XXXK_anal.script for all i1 and i2 directories from mm_300 to mmc_370 but mmc_345".
    # This request: "Do the postprocessing for the results of i1 of each temperature from 295k to 370k". I will assume ALL temps.
    
    temps_int = [295] + list(range(300, 375, 5))
    temps = [f"{t}K" for t in temps_int]
    
    generated_csvs = []
    
    print("Starting Post-Processing for i1 directories...")
    
    for temp in temps:
        xvg_file = get_i1_file(root_dir, temp)
        
        if xvg_file:
            print(f"Processing {temp}: {xvg_file}")
            
            # Run MSD Analysis (creates CSV and D plot)
            try:
                analyze_msd_log.analyze_and_plot(temp, xvg_file)
                
                # Check if CSV was created
                base_name = os.path.splitext(xvg_file)[0]
                csv_path = f"{base_name}.csv"
                if os.path.exists(csv_path):
                    generated_csvs.append(csv_path)
                    
            except Exception as e:
                print(f"Failed analysis for {temp}: {e}")
                
            # Run MSD Plotting (creates raw MSD plot)
            try:
                generate_msd_plots.plot_msd(temp, xvg_file)
            except Exception as e:
                print(f"Failed plotting for {temp}: {e}")
                
        else:
            print(f"Warning: No i1 file found for {temp}")
            
    # Run Arrhenius Plot
    if generated_csvs:
        plot_arrhenius_i1(generated_csvs)
    else:
        print("No CSV data generated to plot Arrhenius.")

if __name__ == "__main__":
    main()
