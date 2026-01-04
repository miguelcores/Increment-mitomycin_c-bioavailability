# Increment Mitomycin C Bioavailability

This project investigates the diffusion properties of Mitomycin C (MMC) in water using Molecular Dynamics (MD) simulations with GROMACS. The goal is to understand the temperature dependence of MMC diffusion to guide research into increasing its bioavailability.

## Simulation Details

*   **System**: 1 molecule of Mitomycin C solvated in a box of 4264 water molecules.
*   **Software**: GROMACS 2021.5.
*   **Temperature Range**: 295K to 370K.

## Repository Structure

The repository is organized by temperature, with specific directories for each simulation run (e.g., `mmc_295K`, `mmc_300K`, etc.).

### Execution Environments

The project supports execution on both High Performance Computing (HPC) clusters and local workstations:

*   **HPC (Slurm)**: Use the `*.script` files located in each temperature directory. These are configured for Slurm workload managers.
*   **Local Workstation (Bash)**: Use the `run_*.sh` scripts (e.g., `run_local.sh`) for running simulations directly in a bash environment without a job scheduler.

## Analysis

Post-processing and analysis are performed using Python scripts included in the root directory.

### Key Scripts

*   **`analyze_msd_log.py`**: Calculates the time-dependent diffusion coefficient $D(t)$ from Mean Squared Displacement (MSD) data. It identifies the diffusive regime, performs log-log analysis, and exports plots and CSV data.
*   **`plot_arrhenius.py`**: Aggregates the calculated diffusion coefficients from all temperatures to generate Arrhenius plots ($\ln(D)$ vs $1/T$) to analyze the temperature dependence of diffusion.
*   **`generate_msd_plots.py`**: Generates standard linear MSD plots.

### Methodology

1.  **MSD Calculation**: GROMACS is used to calculate the MSD from trajectory files.
2.  **Diffusion Coefficient**: $D$ is calculated from the flat region of $D(t) = MSD(t) / (6t)$.
3.  **Arrhenius Analysis**: The relationship between diffusion and temperature is analyzed to check for Arrhenius behavior.
