import os
import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy import stats

# Target files mapping (Temperature -> Filename)
target_files = {
    "295K": "msd_mmc_295K_0-500ns.xvg",
    "300K": "msd_mmc_300K.xvg",
    "305K": "msd_mmc_305K_0-500ns.xvg",
    "310K": "msd_mmc_310K.xvg",
    "315K": "msd_mmc_315K_0-500ns.xvg",
    "320K": "msd_mmc_320K.xvg",
    "325K": "msd_mmc_4265_325K_0-500ns.xvg",
    "330K": "msd_mmc_330K.xvg",
    "335K": "msd_mmc_335K.xvg",
    "340K": "msd_mmc_340K.xvg",
    "345K": "msd_mmc_345K_0-500ns.xvg",
    "350K": "msd_mmc_350K.xvg",
    "355K": "msd_mmc_355K_0-500ns.xvg",
    "360K": "msd_mmc_360K.xvg",
    "365K": "msd_mmc_365K.xvg",
    "370K": "msd_mmc_370K.xvg"
}

def find_csv_files(root_dir, targets):
    found_files = {}
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file.endswith(".csv"):
                # Check if this CSV corresponds to one of our targets
                base_name = os.path.splitext(file)[0]
                # The target values in dict are .xvg files, we need to match base names
                for temp, xvg_name in targets.items():
                    target_base = os.path.splitext(xvg_name)[0]
                    if base_name == target_base:
                        found_files[temp] = os.path.join(root, file)
    return found_files

def parse_csv(filepath):
    data = {}
    with open(filepath, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) == 2:
                data[row[0]] = row[1]
    return data

def parse_d_string(d_str):
    # Format: "0.1234 +/- 0.01" or "N/A"
    if d_str == 'N/A':
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
    csv_files = find_csv_files(root_dir, target_files)
    
    temps = []
    g_ds = []
    g_errs = []
    c_ds = []
    c_errs = []
    
    # Sort by temperature
    sorted_temps = sorted(csv_files.keys(), key=lambda x: int(x.replace('K', '')))
    
    print(f"Found {len(csv_files)} CSV files.")
    
    for temp_str in sorted_temps:
        filepath = csv_files[temp_str]
        data = parse_csv(filepath)
        
        temp_val = float(temp_str.replace('K', ''))
        
        g_d_str = data.get('GROMACS D (1e-5 cm^2/s)', 'N/A')
        c_d_str = data.get('Calculated D (1e-5 cm^2/s)', 'N/A')
        
        g_val, g_err = parse_d_string(g_d_str)
        c_val, c_err = parse_d_string(c_d_str)
        
        if g_val is not None and c_val is not None:
            temps.append(temp_val)
            g_ds.append(g_val)
            g_errs.append(g_err)
            c_ds.append(c_val)
            c_errs.append(c_err)
            print(f"{temp_str}: G_D={g_val}+/-{g_err}, C_D={c_val}+/-{c_err}")
        else:
            print(f"Skipping {temp_str} due to missing data")

    if not temps:
        print("No valid data found.")
        return

    temps = np.array(temps)
    inv_temps = 1000.0 / temps # 1000/T for better scaling
    
    g_ds = np.array(g_ds)
    g_errs = np.array(g_errs)
    c_ds = np.array(c_ds)
    c_errs = np.array(c_errs)
    
    # Log values
    ln_g_ds = np.log(g_ds)
    ln_c_ds = np.log(c_ds)
    
    # Propagated uncertainty: delta(ln x) = delta x / x
    ln_g_errs = g_errs / g_ds
    ln_c_errs = c_errs / c_ds
    
    # Plotting
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: GROMACS D
    ax1.errorbar(inv_temps, ln_g_ds, yerr=ln_g_errs, fmt='o', label='GROMACS D', capsize=5, color='blue')
    slope, intercept, r_value, p_value, std_err = stats.linregress(inv_temps, ln_g_ds)
    fit_line = slope * inv_temps + intercept
    ax1.plot(inv_temps, fit_line, 'r--', label=f'Fit: R²={r_value**2:.4f}\nSlope={slope:.4f}')
    
    ax1.set_title('Arrhenius Plot: GROMACS D')
    ax1.set_xlabel('1000 / T (K⁻¹)')
    ax1.set_ylabel('ln(D)')
    ax1.legend()
    ax1.grid(True)
    
    # Plot 2: Calculated D
    ax2.errorbar(inv_temps, ln_c_ds, yerr=ln_c_errs, fmt='o', label='Calculated D', capsize=5, color='green')
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(inv_temps, ln_c_ds)
    fit_line2 = slope2 * inv_temps + intercept2
    ax2.plot(inv_temps, fit_line2, 'r--', label=f'Fit: R²={r_value2**2:.4f}\nSlope={slope2:.4f}')
    
    ax2.set_title('Arrhenius Plot: Calculated D')
    ax2.set_xlabel('1000 / T (K⁻¹)')
    ax2.set_ylabel('ln(D)')
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig('arrhenius_plots.png')
    print("Saved arrhenius_plots.png")
    plt.close()

if __name__ == "__main__":
    main()
