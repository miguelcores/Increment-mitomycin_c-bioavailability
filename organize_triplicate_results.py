import os
import shutil
import matplotlib.pyplot as plt
import numpy as np

def read_xvg(filepath):
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(('#', '@')):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    data.append([float(parts[0]), float(parts[1])])
                except ValueError:
                    continue
    return np.array(data)

def plot_msd(temp, xvg_path, output_png_path):
    try:
        data = read_xvg(xvg_path)
        if len(data) == 0:
            print(f"Warning: No data in {xvg_path}")
            return

        time_ns = data[:, 0] / 1000.0
        msd_nm2 = data[:, 1]

        plt.figure(figsize=(8, 6))
        plt.plot(time_ns, msd_nm2, label=f'{temp} (Triplicate Avg)', color='blue') # Normal curve
        plt.title(f'Mean Squared Displacement at {temp}\n(Average i0, i1, i2)')
        plt.xlabel('Time (ns)')
        plt.ylabel('MSD (nm²)')
        plt.legend()
        plt.grid(True)
        
        plt.savefig(output_png_path)
        plt.close()
        print(f"Generated plot: {output_png_path}")
    except Exception as e:
        print(f"Error plotting {xvg_path}: {e}")

def main():
    root_dir = os.getcwd()
    base_folder_name = "average_results_i0_i1_i2"
    base_folder_path = os.path.join(root_dir, base_folder_name)
    
    if not os.path.exists(base_folder_path):
        print(f"Error: Folder '{base_folder_name}' not found.")
        return

    print(f"Organizing results in {base_folder_path}...")

    temps_int = list(range(300, 375, 5))
    
    for t in temps_int:
        temp_str = f"{t}K"
        
        # Define the subdirectory for this temperature: average_results_i0_i1_i2/mmc_300K
        temp_subdir_name = f"mmc_{temp_str}"
        temp_subdir_path = os.path.join(base_folder_path, temp_subdir_name)
        
        if not os.path.exists(temp_subdir_path):
            os.makedirs(temp_subdir_path)
            print(f"Created directory: {temp_subdir_path}")

        # Files to move: xvg, csv, and existing plots
        files_to_move = [
            f"msd_mmc_{temp_str}.xvg", 
            f"msd_mmc_{temp_str}.csv",
            f"msd_mmc_{temp_str}_D_vs_t.png",
            f"msd_mmc_{temp_str}_log.png"
        ]
        
        xvg_file_path = None

        for fname in files_to_move:
            src_path = os.path.join(base_folder_path, fname)
            dst_path = os.path.join(temp_subdir_path, fname)
            
            # Check if file is in root of base_folder
            if os.path.exists(src_path):
                shutil.move(src_path, dst_path)
                print(f"Moved {fname} -> {temp_subdir_name}/")
                if fname.endswith(".xvg"):
                    xvg_file_path = dst_path
            # Check if file is already in destination (e.g. from previous run)
            elif os.path.exists(dst_path):
                if fname.endswith(".xvg"):
                    xvg_file_path = dst_path
        
        # Generate plot if we have the xvg file
        if xvg_file_path and os.path.exists(xvg_file_path):
            plot_filename = f"msd_triplicate_avg_plot_{temp_str}.png"
            plot_path = os.path.join(temp_subdir_path, plot_filename)
            plot_msd(temp_str, xvg_file_path, plot_path)
        else:
            print(f"Skipping plot for {temp_str}: XVG file not found.")

    print("Done organizing and plotting.")

if __name__ == "__main__":
    main()
