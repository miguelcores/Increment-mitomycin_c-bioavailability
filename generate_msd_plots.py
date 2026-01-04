import os
import matplotlib.pyplot as plt
import numpy as np

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

def read_xvg(filepath):
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(('#', '@')):
                continue
            parts = line.split()
            # Assuming 1st column is time (ps), 2nd is MSD
            if len(parts) >= 2:
                try:
                    data.append([float(parts[0]), float(parts[1])])
                except ValueError:
                    continue
    return np.array(data)

def find_files(root_dir, targets):
    found_files = {}
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file in targets.values():
                # Check if we already found this file (to avoid duplicates if any, or handle multiple locations)
                # For now, let's assume unique filenames or take the first one found
                # We need to map back to temperature
                for temp, target_name in targets.items():
                    if file == target_name:
                        found_files[temp] = os.path.join(root, file)
    return found_files

def plot_msd(temp, filepath):
    try:
        data = read_xvg(filepath)
        if len(data) == 0:
            print(f"No data found in {filepath}")
            return

        time = data[:, 0] / 1000.0 # Convert ps to ns
        msd = data[:, 1]

        plt.figure(figsize=(8, 6))
        plt.plot(time, msd, label=f'{temp}')
        plt.title(f'Mean Squared Displacement at {temp}')
        plt.xlabel('Time (ns)')
        plt.ylabel('MSD (nm²)')
        plt.legend()
        plt.grid(True)

        output_path = os.path.splitext(filepath)[0] + '.png'
        plt.savefig(output_path)
        plt.close()
        print(f"Saved plot to {output_path}")

    except Exception as e:
        print(f"Error processing {filepath}: {e}")

def main():
    root_dir = os.getcwd()
    print(f"Searching for files in {root_dir}...")
    
    found_files = find_files(root_dir, target_files)
    
    for temp, filename in target_files.items():
        if temp in found_files:
            print(f"Processing {temp}: {found_files[temp]}")
            plot_msd(temp, found_files[temp])
        else:
            print(f"Warning: File for {temp} ({filename}) not found.")

if __name__ == "__main__":
    main()
