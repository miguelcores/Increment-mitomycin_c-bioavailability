# Mitomycin C Bioavailability: Computational Chemistry & Drug Delivery Research

A comprehensive computational and experimental study of mitomycin C (MMC) molecular behavior to enhance drug bioavailability through temperature-dependent kinetics, concentration-dependent effects, and membrane interaction analysis.

## 🎯 Project Overview

This research combines:
- **Molecular Dynamics Simulations**: Temperature and concentration studies of MMC diffusion
- **Membrane Biophysics**: DPPC lipid bilayer interactions with MMC
- **Life Cycle Assessment**: Environmental impact analysis of drug formulation
- **Collaborative Analysis**: Integration of computational and experimental data

## 📊 Research Components

### 1. Temperature Study (295K-370K)
Analysis of MMC diffusion kinetics across physiologically relevant temperatures.

**Location**: `temperature_study/`

- Multiple systems at different temperatures
- Measure Mean Square Displacement (MSD) to extract diffusion coefficients
- Generate Arrhenius plots to determine activation energy
- Expected outcomes: Temperature-dependent diffusion parameters

**Key Files**:
- `*/run_*_mmc.script` - SLURM submission scripts for MD simulations
- `*/run_*_mmc_anal.script` - Analysis job scripts for trajectory processing

### 2. Concentration Study (100-180 MMC molecules)
Investigation of concentration-dependent MMC behavior at T=325K.

**Location**: `concentration_study/`

- `300K/` - 300 Kelvin concentration screening
- `325K/` - 325 Kelvin comprehensive concentration series
- Systems from 2 to 152+ MMC molecules per simulation box

**Key Files**:
- `*/run_*_mmc.script` - Production MD simulations
- `*/run_*_mmc_anal.script` - Trajectory analysis scripts
- `*.csv` - Analyzed density vs concentration results

### 3. DPPC Membrane System Studies
Characterization of MMC interactions with dipalmitoylphosphatidylcholine (DPPC) lipid bilayers.

**Location**: `DPPC/`

- Multiple replicas (i0, i1) for each system
- 1 to 100+ MMC molecules interacting with DPPC membranes
- Parallel replicas for statistical robustness

**Key Files**:
- `topol.top` - System topology files (modified for MMC-DPPC)
- `mdp/` - GROMACS molecular dynamics parameters
- `*/i0/`, `*/i1/` - Replica directories with simulation data

### 4. Mitomycin C Force Field
Custom force field parameters for MMC in GROMOS force field.

**Location**: `mitomycin_c/`

- `mmc.gro` - Mitomycin C structure file
- `mmc.itp` - Force field parameters (ATB-derived)
- `gromos54a7_atb.ff/` - Complete force field library

### 5. Environmental Impact Analysis
Life cycle assessment and carbon footprint of formulation strategies.

**Location**: `arch/` (archived documentation)

- `ACV_Ambiental_Mitomicina_C.ipynb` - Environmental analysis notebook
- Carbon footprint calculations

## 🔬 Methodology

### Simulation Details

**Software Stack**:
- GROMACS 2021.5 - Molecular Dynamics engine
- GROMOS54a7 force field - for MMC and water
- OpenMPI - Parallel computing framework

**Computational Resources**:
- HPC cluster submission via SLURM
- Typical runs: 48-96 hour walltime
- 1-48 CPU cores per simulation

### Analysis Workflows

**For Temperature Studies**:
1. Run MD simulations at each temperature
2. Extract center-of-mass coordinates
3. Calculate MSD as function of time
4. Fit diffusion coefficients (Einstein relation)
5. Generate Arrhenius plots
6. Extract activation energy (Ea)

**For Concentration Studies**:
1. Build systems with increasing MMC count
2. Run equilibration and production phases
3. Analyze density fluctuations
4. Detect aggregation/phase transitions
5. Quantify concentration effects on solubility

**For Membrane Studies**:
1. Insert MMC into DPPC bilayer
2. Monitor adsorption/desorption kinetics
3. Measure membrane perturbation
4. Calculate interaction energies
5. Determine permeability rates

## 📁 Directory Structure

```
.
├── README.md                           # This file
├── concentration_study/
│   ├── 300K/                          # 300K concentration screening
│   ├── 325K/                          # 325K main concentration series
│   ├── analyze_density_vs_concentration.py
│   └── hpc_*.sh                       # HPC submission scripts
├── DPPC/
│   ├── [N]_mmc/                       # N = 1,2,4,8,16,32,50,64,100
│   │   ├── i0/, i1/                   # Replica 0 and 1
│   │   ├── topol.top
│   │   └── run_*.script
│   ├── mdp/                           # MD parameter files
│   ├── gromos53a6_lipid.ff/           # DPPC force field
│   └── topol.top
├── temperature_study/
│   ├── [N]_mmc/
│   │   ├── [TEMP]K/
│   │   │   ├── i0/, i1/
│   │   │   ├── run_*_mmc.script
│   │   │   └── run_*_mmc_anal.script
│   └── submit_all_[TEMP]K_slurm.sh
├── mitomycin_c/
│   ├── mmc.gro
│   ├── mmc.itp
│   └── gromos54a7_atb.ff/
├── arch/                              # Project documentation
│   ├── QUICK_START_48_HOURS.md
│   ├── RESEARCH_STRATEGY_AND_PAPER_STRUCTURE.md
│   ├── SIMULATION_QUALITY_DIAGNOSTICS.md
│   ├── COLLABORATION_GUIDE.md
│   └── ACV_Ambiental_Mitomicina_C.ipynb
├── carbon_footprint.txt               # Environmental impact summary
└── useful_commands.txt                # Archived command reference
```

## 🚀 Quick Start

### Prerequisites
- GROMACS 2021.5+
- OpenMPI for parallel runs
- Python 3.8+ for analysis scripts
- HPC cluster access with SLURM scheduler

### Running a Simulation

1. **Navigate to target system**:
   ```bash
   cd concentration_study/325K/32_mmc
   ```

2. **Submit the job**:
   ```bash
   sbatch run_32_mmc.script
   ```

3. **Monitor progress**:
   ```bash
   squeue -u $USER
   tail -f slurm-*.out
   ```

4. **Analyze results**:
   ```bash
   sbatch run_32_mmc_anal.script
   ```

### Post-Processing Analysis

See `arch/QUICK_START_48_HOURS.md` for detailed analysis workflows including:
- Arrhenius plot generation
- MSD calculations
- Diffusion coefficient extraction
- Quality diagnostics

## 📈 Expected Outputs

### From Temperature Studies
- Diffusion coefficients (D) at 295K-370K
- Arrhenius plots with activation energy (Ea)
- Temperature-dependent transport properties

### From Concentration Studies
- Density-concentration relationships
- Critical concentration for aggregation
- Solubility limits
- Phase behavior analysis

### From Membrane Studies
- MMC permeation rates across DPPC
- Adsorption/desorption kinetics
- Membrane perturbation parameters
- Cellular uptake insights

## 📚 Documentation

- **`QUICK_START_48_HOURS.md`** - Rapid setup for newcomers
- **`RESEARCH_STRATEGY_AND_PAPER_STRUCTURE.md`** - Research roadmap and publication strategy
- **`SIMULATION_QUALITY_DIAGNOSTICS.md`** - Troubleshooting and data quality checks
- **`COLLABORATION_GUIDE.md`** - Instructions for inter-group collaboration
- **`ACV_Ambiental_Mitomicina_C.ipynb`** - Environmental assessment notebook

## 📋 Key Analysis Scripts

- `create_slurm_scripts.py` - Generate SLURM submission scripts
- `concentration_study/*/run_*_mmc.script` - MD production runs
- `concentration_study/*/run_*_mmc_anal.script` - Post-simulation analysis
- `DPPC/build_dppc_mmc_hydrated_cases.sh` - System preparation for membrane studies

## 🔗 Force Field References

- **GROMOS54a7**: For protein/drug parameterization
- **GROMOS53a6 Lipid**: For DPPC membrane
- **ATB (Automated Topology Builder)**: Source for MMC parameters

## 📊 Citation & Attribution

This project integrates:
1. Computational chemistry (MMC diffusion kinetics)
2. Membrane biophysics (DPPC interactions)
3. Environmental impact assessment
4. Collaborative experimental data

For publications, please cite:
- GROMACS: Abraham et al. (2015)
- Force fields: Original publications of GROMOS
- This project: [Your Publication Details]

## 🤝 Collaboration

This work involves collaborative efforts between:
- Computational chemistry team (MD simulations)
- Biophysics group (membrane experiments)
- Environmental science team (LCA analysis)

See `arch/COLLABORATION_GUIDE.md` for coordination details.

## ⚠️ Important Notes

- **Large Data Files**: Trajectory files (*.trr, *.xtc) are excluded from version control
- **Sensitive Information**: Personal credentials removed; use local HPC configuration
- **Environment**: Assumes Linux/Unix environment with SLURM; WSL for Windows users
- **Dependencies**: See module loads in submission scripts for required software

## 📝 License

[Specify your license here - e.g., MIT, GPL, Creative Commons, etc.]

## ✉️ Questions?

For questions about methodology, analysis, or collaboration:
- Review `arch/QUICK_START_48_HOURS.md`
- Check `SIMULATION_QUALITY_DIAGNOSTICS.md` for troubleshooting
- Consult `arch/COLLABORATION_GUIDE.md` for coordination

---

**Last Updated**: May 2026  
**Status**: Active Research  
**Version**: Final Release
