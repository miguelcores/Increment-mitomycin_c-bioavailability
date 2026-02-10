import os
import matplotlib.pyplot as plt
import numpy as np
import re
import csv

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

def parse_xvg(filepath):
    data = []
    gromacs_d = None
    gromacs_err = None
    
    # Regex to find D value in comments
    # Example: # D[     MMC] = 0.1234 (+/- 0.01) (1e-5 cm^2/s)
    d_pattern = re.compile(r"D\[.*\]\s*=\s*([0-9\.]+)\s*\(\+/-\s*([0-9\.]+)\)")

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Check for GROMACS D in header
            if line.startswith("#") and "D[" in line:
                match = d_pattern.search(line)
                if match:
                    try:
                        gromacs_d = float(match.group(1))
                        gromacs_err = float(match.group(2))
                    except ValueError:
                        pass

            if not line or line.startswith(('#', '@')):
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                try:
                    data.append([float(parts[0]), float(parts[1])])
                except ValueError:
                    continue
    return np.array(data), gromacs_d, gromacs_err

def find_files(root_dir, targets):
    found_files = {}
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file in targets.values():
                for temp, target_name in targets.items():
                    if file == target_name:
                        found_files[temp] = os.path.join(root, file)
    return found_files

def analyze_and_plot(temp, filepath):
    try:
        data, gromacs_d, gromacs_err = parse_xvg(filepath)
        if len(data) == 0:
            print(f"No data found in {filepath}")
            return

        time = data[:, 0] # ps
        msd = data[:, 1]  # nm^2

        # Filter out zero or negative values
        valid_indices = (time > 0) & (msd > 0)
        time = time[valid_indices]
        msd = msd[valid_indices]
        
        # Calculate D(t) = MSD(t) / (6 * t)
        # Units: nm^2 / ps
        d_t = msd / (6.0 * time)
        
        # Find the "flat" region of D(t)
        # Heuristic:
        # 1. Ignore the first 10% of the trajectory (equilibration/ballistic)
        # 2. Use a sliding window of 30% of the total length
        # 3. Find the window with the minimum relative standard deviation (CV)
        
        n_points = len(time)
        start_idx_search = int(n_points * 0.1)
        window_size = int(n_points * 0.3)
        
        if window_size < 10:
            window_size = n_points // 2 # Fallback for very short files
            
        best_cv = float('inf')
        best_start_idx = start_idx_search
        best_end_idx = start_idx_search + window_size
        
        # Sliding window
        # We step by 1% of points to speed up
        step = max(1, int(n_points * 0.01))
        
        for i in range(start_idx_search, n_points - window_size, step):
            segment = d_t[i : i + window_size]
            mean_val = np.mean(segment)
            std_val = np.std(segment)
            
            if mean_val > 0:
                cv = std_val / mean_val
                if cv < best_cv:
                    best_cv = cv
                    best_start_idx = i
                    best_end_idx = i + window_size

        # Calculate final D from the best window
        d_raw_segment = d_t[best_start_idx : best_end_idx]
        d_raw = np.mean(d_raw_segment) # nm^2/ps
        d_raw_std = np.std(d_raw_segment) # nm^2/ps
        
        d_calc = d_raw * 1000.0 # Convert to 10^-5 cm^2/s
        d_calc_err = d_raw_std * 1000.0 # Convert to 10^-5 cm^2/s
        
        # --- Plot 1: D(t) vs t ---
        plt.figure(figsize=(10, 6))
        plt.plot(time, d_t * 1000.0, label='D(t) = MSD/(6t)', color='blue', alpha=0.6)
        
        # Highlight the selected region
        fit_time = time[best_start_idx : best_end_idx]
        fit_d = d_t[best_start_idx : best_end_idx] * 1000.0
        plt.plot(fit_time, fit_d, color='red', linewidth=2, label='Selected Diffusive Region')
        plt.axhline(y=d_calc, color='green', linestyle='--', label=f'Avg D = {d_calc:.4f} +/- {d_calc_err:.4f}')
        
        plt.title(f'Diffusion Coefficient D(t) vs Time at {temp}\nCalculated D = {d_calc:.4f} +/- {d_calc_err:.4f} (1e-5 cm^2/s)')
        plt.xlabel('Time (ps)')
        plt.ylabel('D (1e-5 cm^2/s)')
        plt.legend()
        plt.grid(True)
        
        output_dir = os.path.dirname(filepath)
        base_name = os.path.splitext(os.path.basename(filepath))[0]
        
        plot_d_filename = f"{base_name}_D_vs_t.png"
        plot_d_path = os.path.join(output_dir, plot_d_filename)
        plt.savefig(plot_d_path)
        plt.close()
        print(f"Saved D(t) plot to {plot_d_path}")

        # --- Plot 2: Log-Log MSD vs t ---
        # Calculate actual slope (alpha) in the selected region
        fit_time = time[best_start_idx : best_end_idx]
        fit_msd_data = msd[best_start_idx : best_end_idx]
        
        log_fit_time = np.log(fit_time)
        log_fit_msd = np.log(fit_msd_data)
        
        slope, intercept = np.polyfit(log_fit_time, log_fit_msd, 1)
        alpha = slope
        
        plt.figure(figsize=(10, 8))
        plt.loglog(time, msd, '.', label='Data', markersize=2, alpha=0.5)
        
        # Plot the actual fit line
        fit_msd_actual = np.exp(intercept) * (fit_time ** slope)
        
        label_text = (f'Fit Region: {fit_time[0]:.0f}-{fit_time[-1]:.0f} ps\n'
                      f'Slope (α) = {alpha:.4f}')
        
        plt.loglog(fit_time, fit_msd_actual, 'r-', linewidth=2, label=label_text)
        
        plt.title(f'MSD Log-Log Plot at {temp}')
        plt.xlabel('Time (ps)')
        plt.ylabel('MSD (nm^2)')
        plt.legend()
        plt.grid(True, which="both", ls="-", alpha=0.5)
        
        plot_log_filename = f"{base_name}_log.png"
        plot_log_path = os.path.join(output_dir, plot_log_filename)
        plt.savefig(plot_log_path)
        plt.close()
        print(f"Saved log-log plot to {plot_log_path}")
        
        # Save CSV
        csv_filename = f"{base_name}.csv"
        csv_path = os.path.join(output_dir, csv_filename)
        
        with open(csv_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Parameter', 'Value'])
            writer.writerow(['Temperature', temp])
            
            g_d_str = f"{gromacs_d} +/- {gromacs_err}" if gromacs_d is not None else 'N/A'
            writer.writerow(['GROMACS D (1e-5 cm^2/s)', g_d_str])
            
            c_d_str = f"{d_calc:.4f} +/- {d_calc_err:.4f}"
            writer.writerow(['Calculated D (1e-5 cm^2/s)', c_d_str])
            
            writer.writerow(['Method', 'Flat region of MSD/(6t)'])
            writer.writerow(['Fit Start Time (ps)', time[best_start_idx]])
            writer.writerow(['Fit End Time (ps)', time[best_end_idx]])
            writer.writerow(['Fitted Slope (alpha)', alpha])
            
        print(f"Saved CSV to {csv_path}")

    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        import traceback
        traceback.print_exc()

def main():
    root_dir = os.getcwd()
    print(f"Searching for files in {root_dir}...")
    
    found_files = find_files(root_dir, target_files)
    
    for temp, filename in target_files.items():
        if temp in found_files:
            print(f"Processing {temp}: {found_files[temp]}")
            analyze_and_plot(temp, found_files[temp])
        else:
            print(f"Warning: File for {temp} ({filename}) not found.")

if __name__ == "__main__":
    main()
