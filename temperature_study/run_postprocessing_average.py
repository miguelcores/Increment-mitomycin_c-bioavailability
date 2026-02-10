import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import csv
import re

# Import existing modules
import analyze_msd_log
import run_postprocessing_i1
import run_postprocessing_i2

def parse_xvg_header_d(filepath):
    """
    Extracts GROMACS D from the file header directly.
    """
    d_pattern = re.compile(r"D\[.*\]\s*=\s*([0-9\.]+)\s*\(\+/-\s*([0-9\.]+)\)")
    val, err = None, None
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("#") and "D[" in line:
                match = d_pattern.search(line)
                if match:
                    try:
                        val = float(match.group(1))
                        err = float(match.group(2))
                    except ValueError:
                        pass
                break # Found it
    return val, err

def create_averaged_xvg(temp, file1, file2, output_path):
    # Parse data using analyze_msd_log's function
    data1, d1, err1 = analyze_msd_log.parse_xvg(file1)
    data2, d2, err2 = analyze_msd_log.parse_xvg(file2)
    
    if len(data1) == 0 or len(data2) == 0:
        print(f"Error: Empty data for {temp}")
        return False
        
    # Check lengths
    min_len = min(len(data1), len(data2))
    data1 = data1[:min_len]
    data2 = data2[:min_len]
    
    # Check times align
    if not np.allclose(data1[:, 0], data2[:, 0], atol=1e-3):
        print(f"Error: Time mismatch for {temp}")
        return False
        
    # Average MSD
    avg_data = np.copy(data1)
    avg_data[:, 1] = (data1[:, 1] + data2[:, 1]) / 2.0
    
    # Average GROMACS D calculations if present
    avg_d = None
    avg_err = None
    if d1 is not None and d2 is not None:
        avg_d = (d1 + d2) / 2.0
        # Error propagation for average: sqrt(err1^2 + err2^2) / 2
        avg_err = np.sqrt(err1**2 + err2**2) / 2.0
        
    # Write new XVG file
    with open(output_path, 'w') as f:
        f.write("# Averaged MSD from i1 and i2\n")
        f.write(f"# Created by run_postprocessing_average.py\n")
        if avg_d is not None:
             f.write(f"# D[ MMC] = {avg_d:.4f} (+/- {avg_err:.4f}) (1e-5 cm^2/s)\n")
        f.write("@    title \"Averaged Mean Square Displacement\"\n")
        f.write("@    xaxis  label \"Time (ps)\"\n")
        f.write("@    yaxis  label \"MSD (nm^2)\"\n")
        f.write("@TYPE xy\n")
        
        for row in avg_data:
            f.write(f"{row[0]:.4f}  {row[1]:.6f}\n")
            
    return True

def plot_arrhenius_avg(csv_files, output_filename='arrhenius_plots_average.png'):
    temps = []
    g_ds = []
    g_errs = []
    c_ds = []
    c_errs = []
    
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

    print("\nArrhenius Data (Averaged i1/i2):")
    for d in data_list:
        print(f"{d['T']}K: Calc D = {d['C_D']:.4f}, Gromacs Avg D = {d['G_D']}")
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
        ax1.errorbar(g_inv_temps, ln_g_ds, yerr=ln_g_errs, fmt='o', label='Avg GROMACS D', capsize=5, color='blue')
        slope, intercept, r_value, _, _ = stats.linregress(g_inv_temps, ln_g_ds)
        fit_line = slope * g_inv_temps + intercept
        ax1.plot(g_inv_temps, fit_line, 'r--', label=f'Fit: R²={r_value**2:.4f}\nSlope={slope:.4f}')
        ax1.set_title('Arrhenius: Averaged GROMACS D')
        ax1.set_xlabel('1000 / T (K⁻¹)')
        ax1.set_ylabel('ln(D)')
        ax1.legend()
        ax1.grid(True)
    
    # Plot 2: Calculated Avg D
    ax2.errorbar(c_inv_temps, ln_c_ds, yerr=ln_c_errs, fmt='o', label='Avg Calculated D', capsize=5, color='green')
    slope, intercept, r_value, _, _ = stats.linregress(c_inv_temps, ln_c_ds)
    fit_line = slope * c_inv_temps + intercept
    ax2.plot(c_inv_temps, fit_line, 'r--', label=f'Fit: R²={r_value**2:.4f}\nSlope={slope:.4f}')
    ax2.set_title('Arrhenius: Averaged Calculated D')
    ax2.set_xlabel('1000 / T (K⁻¹)')
    ax2.set_ylabel('ln(D)')
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig(output_filename)
    plt.close()
    print(f"\nSaved averaged Arrhenius plot to {output_filename}")


def main():
    root_dir = os.getcwd()
    avg_dir = os.path.join(root_dir, "average_results")
    if not os.path.exists(avg_dir):
        os.makedirs(avg_dir)
        
    temps_int = list(range(300, 375, 5))
    temps = [f"{t}K" for t in temps_int]
    
    generated_csvs = []
    
    print("Starting Averaged Post-Processing (i1 + i2)...")
    
    for temp in temps:
        file1 = run_postprocessing_i1.get_i1_file(root_dir, temp)
        file2 = run_postprocessing_i2.get_i2_file(root_dir, temp)
        
        if file1 and file2:
            print(f"Processing {temp}...")
            # Create averaged file
            avg_xvg_path = os.path.join(avg_dir, f"msd_mmc_{temp}.xvg")
            success = create_averaged_xvg(temp, file1, file2, avg_xvg_path)
            
            if success:
                # Run standard analysis on the averaged file
                analyze_msd_log.analyze_and_plot(temp, avg_xvg_path)
                
                # Check for CSV
                csv_path = os.path.splitext(avg_xvg_path)[0] + ".csv"
                if os.path.exists(csv_path):
                    generated_csvs.append(csv_path)
        else:
            print(f"Skipping {temp}: Missing i1 or i2 file")
            if not file1: print(f"  Missing i1")
            if not file2: print(f"  Missing i2")
            
    # Arrhenius
    if generated_csvs:
        plot_arrhenius_avg(generated_csvs)


if __name__ == "__main__":
    main()
