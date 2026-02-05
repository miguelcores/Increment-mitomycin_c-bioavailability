import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import csv
import re
import glob

# Import existing modules
import analyze_msd_log
import run_postprocessing_i1
import run_postprocessing_i2

def get_original_file(root_dir, temp):
    """
    Finds the 'original' / 'first postprocessing' file.
    Priority: i0 > origin > root of temp folder.
    """
    target_name = analyze_msd_log.target_files.get(temp)
    
    # Locate temp directory
    temp_num = temp.replace('K', '')
    pattern = os.path.join(root_dir, f"mmc_{temp_num}[Kk]*")
    dirs = glob.glob(pattern)
    if not dirs: return None
    
    # Use the first matching directory (e.g. mmc_300K)
    temp_dir = dirs[0]
    
    # Priority locations
    subdirs_to_check = ['i0', 'origin', '']
    
    for sub in subdirs_to_check:
        search_dir = os.path.join(temp_dir, sub)
        if not os.path.isdir(search_dir):
            continue
            
        # 1. Try exact target name from analyze_msd_log dict
        if target_name:
            cand = os.path.join(search_dir, target_name)
            if os.path.exists(cand):
                return cand
        
        # 2. Try generic pattern "msd_mmc_*.xvg"
        pattern = os.path.join(search_dir, "msd_mmc_*.xvg")
        cands = glob.glob(pattern)
        if cands:
             # Just take the first valid xvg found
             # Since glob doesn't recurse by default, we won't accidentally pick up i1/i2 files if we are scanning root
             return cands[0]
             
    return None

def create_triplicate_averaged_xvg(temp, files, output_path):
    """
    Averages MSD columns from a list of files.
    """
    parsed_data = [] # List of (data_array, d_val, err_val)
    
    for fpath in files:
        if fpath and os.path.exists(fpath):
            d_arr, d_val, d_err = analyze_msd_log.parse_xvg(fpath)
            if len(d_arr) > 0:
                parsed_data.append((d_arr, d_val, d_err))
            else:
                print(f"Warning: Empty data in {fpath}")
        else:
             print(f"Warning: File missing: {fpath}")

    if not parsed_data:
        print(f"Error: No valid data for {temp}")
        return False
        
    # Determine minimum length to truncate to
    min_len = min(len(x[0]) for x in parsed_data)
    
    # Valid datasets aligned by time?
    # Check first file's time vs others
    ref_time = parsed_data[0][0][:min_len, 0]
    
    for i in range(1, len(parsed_data)):
        curr_time = parsed_data[i][0][:min_len, 0]
        if not np.allclose(ref_time, curr_time, atol=1e-3):
            print(f"Error: Time mismatch for {temp} in file {i+1}")
            # Try to align? For now, strict fail or just truncate mismatch
            # We'll just return False to be safe
            return False

    # Average MSD
    # Stack MSD columns: shape (min_len, n_files)
    msd_stack = np.column_stack([d[0][:min_len, 1] for d in parsed_data])
    avg_msd = np.mean(msd_stack, axis=1)
    
    # Construct averaged data array (Time, MSD)
    avg_data = np.column_stack((ref_time, avg_msd))
    
    # Average GROMACS D if available
    # Collect available Ds
    valid_ds = [d[1] for d in parsed_data if d[1] is not None]
    valid_errs = [d[2] for d in parsed_data if d[2] is not None]
    
    avg_d = None
    avg_err = None
    if valid_ds:
        avg_d = np.mean(valid_ds)
        # Propagate error: sqrt(sum(err^2))/N
        if valid_errs and len(valid_errs) == len(valid_ds):
            avg_err = np.sqrt(np.sum(np.array(valid_errs)**2)) / len(valid_ds)
            
    # Write new XVG file
    with open(output_path, 'w') as f:
        f.write("# Averaged MSD from i0/Original, i1, and i2\n")
        f.write(f"# Created by run_postprocessing_average_triplicate.py\n")
        f.write(f"# Averaged over {len(parsed_data)} replicates\n")
        if avg_d is not None:
             f.write(f"# D[ MMC] = {avg_d:.4f} (+/- {avg_err:.4f}) (1e-5 cm^2/s)\n")
        f.write("@    title \"Averaged Mean Square Displacement (Triplicate)\"\n")
        f.write("@    xaxis  label \"Time (ps)\"\n")
        f.write("@    yaxis  label \"MSD (nm^2)\"\n")
        f.write("@TYPE xy\n")
        
        for row in avg_data:
            f.write(f"{row[0]:.4f}  {row[1]:.6f}\n")
            
    return True

def plot_arrhenius_triplicate(csv_files, output_filename='arrhenius_plots_triplicate_avg.png'):
    temps = []
    data_list = []
    
    # Read CSVs
    for filepath in csv_files:
        try:
            with open(filepath, 'r') as f:
                reader = csv.reader(f)
                file_data = {}
                for row in reader:
                    if len(row) == 2:
                        file_data[row[0]] = row[1]
                
                temp_str = file_data.get('Temperature', 'N/A')
                # Clean temp string
                try:
                    t_val = float(re.sub(r'[^\d\.]', '', temp_str))
                except:
                    continue
                
                g_d_val, g_d_err = run_postprocessing_i1.parse_d_string(file_data.get('GROMACS D (1e-5 cm^2/s)', 'N/A'))
                c_d_val, c_d_err = run_postprocessing_i1.parse_d_string(file_data.get('Calculated D (1e-5 cm^2/s)', 'N/A'))
                
                if c_d_val is not None:
                    data_list.append({
                        'T': t_val,
                        'G_D': g_d_val, 'G_Err': g_d_err,
                        'C_D': c_d_val, 'C_Err': c_d_err
                    })
        except Exception:
            pass

    data_list.sort(key=lambda x: x['T'])
    
    calc_data = []
    grom_data = []

    print("\nArrhenius Data (Triplicate Average):")
    for d in data_list:
        val_str = f"{d['C_D']:.4f}"
        grom_str = f"{d['G_D']:.4f}" if d['G_D'] is not None else "None"
        print(f"{d['T']}K: Calc D = {val_str}, Gromacs Avg D = {grom_str}")
        
        calc_data.append((d['T'], d['C_D'], d['C_Err']))
        if d['G_D'] is not None:
             grom_data.append((d['T'], d['G_D'], d['G_Err']))

    if not calc_data:
        print("No valid data for plotting.")
        return

    # Prepare arrays
    c_temps = np.array([x[0] for x in calc_data])
    c_inv_temps = 1000.0 / c_temps
    c_ds = np.array([x[1] for x in calc_data])
    c_errs = np.array([x[2] for x in calc_data])
    ln_c_ds = np.log(c_ds)
    ln_c_errs = c_errs / c_ds

    if grom_data:
        g_temps = np.array([x[0] for x in grom_data])
        g_inv_temps = 1000.0 / g_temps
        g_ds = np.array([x[1] for x in grom_data])
        g_errs = np.array([x[2] for x in grom_data])
        ln_g_ds = np.log(g_ds)
        ln_g_errs = g_errs / g_ds
    else:
        g_inv_temps = None

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: GROMACS Avg D
    if grom_data:
        ax1.errorbar(g_inv_temps, ln_g_ds, yerr=ln_g_errs, fmt='o', label='Tri-Avg GROMACS D', capsize=5, color='blue')
        slope, intercept, r_value, _, _ = stats.linregress(g_inv_temps, ln_g_ds)
        fit_line = slope * g_inv_temps + intercept
        ax1.plot(g_inv_temps, fit_line, 'r--', label=f'Fit: R²={r_value**2:.4f}\nSlope={slope:.4f}')
        ax1.set_title('Arrhenius: Tri-Averaged GROMACS D')
        ax1.set_xlabel('1000 / T (K⁻¹)')
        ax1.set_ylabel('ln(D)')
        ax1.legend()
        ax1.grid(True)
    else:
        ax1.text(0.5, 0.5, "No Data", ha='center')
    
    # Plot 2: Calculated Avg D
    ax2.errorbar(c_inv_temps, ln_c_ds, yerr=ln_c_errs, fmt='o', label='Tri-Avg Calculated D', capsize=5, color='green')
    slope, intercept, r_value, _, _ = stats.linregress(c_inv_temps, ln_c_ds)
    fit_line = slope * c_inv_temps + intercept
    ax2.plot(c_inv_temps, fit_line, 'r--', label=f'Fit: R²={r_value**2:.4f}\nSlope={slope:.4f}')
    ax2.set_title('Arrhenius: Tri-Averaged Calculated D')
    ax2.set_xlabel('1000 / T (K⁻¹)')
    ax2.set_ylabel('ln(D)')
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig(output_filename)
    plt.close()
    print(f"\nSaved triplicate averaged Arrhenius plot to {output_filename}")

def main():
    root_dir = os.getcwd()
    avg_dir = os.path.join(root_dir, "average_triplicate_results")
    if not os.path.exists(avg_dir):
        os.makedirs(avg_dir)
        
    # User requested 300K to 370K
    temps_int = list(range(300, 375, 5))
    temps = [f"{t}K" for t in temps_int]
    
    generated_csvs = []
    
    print("Starting Post-Processing for Triplicate Average (Original + i1 + i2)...")
    
    for temp in temps:
        # Find files
        f0 = get_original_file(root_dir, temp)
        f1 = run_postprocessing_i1.get_i1_file(root_dir, temp)
        f2 = run_postprocessing_i2.get_i2_file(root_dir, temp)
        
        # We process if at least one file exists, though quality improves with more
        files = [f for f in [f0, f1, f2] if f is not None]
        
        print(f"Processing {temp}: Found {len(files)} files.")
        if f0: print(f"  - Original: {os.path.basename(f0)}")
        else: print(f"  - Original: [MISSING]")
        
        if f1: print(f"  - i1: {os.path.basename(f1)}")
        else: print(f"  - i1: [MISSING]")

        if f2: print(f"  - i2: {os.path.basename(f2)}")
        else: print(f"  - i2: [MISSING]")
        
        if len(files) > 0:
            # Create averaged file
            avg_xvg_path = os.path.join(avg_dir, f"msd_mmc_{temp}.xvg")
            success = create_triplicate_averaged_xvg(temp, files, avg_xvg_path)
            
            if success:
                # Run standard analysis
                analyze_msd_log.analyze_and_plot(temp, avg_xvg_path)
                
                # Check for CSV
                csv_path = os.path.splitext(avg_xvg_path)[0] + ".csv"
                if os.path.exists(csv_path):
                    generated_csvs.append(csv_path)
        else:
            print(f"Skipping {temp}: No files found.")

    # Arrhenius
    if generated_csvs:
        plot_arrhenius_triplicate(generated_csvs)

if __name__ == "__main__":
    main()
