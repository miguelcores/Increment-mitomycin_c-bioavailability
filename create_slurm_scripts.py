import os

def main():
    root_dir = os.getcwd()
    conc_study_dir = os.path.join(root_dir, "concentration_study")
    
    concentrations = [2, 4, 8, 12, 16]
    temperature = "325K"
    
    template = """#!/bin/bash

##----------------------- Start job description -----------------------

#SBATCH --partition=standard
#SBATCH --job-name=mmc_{n_mmc}_325K
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@institution.org

##----------------------- End job description ------------------------

module load apps/2021
module load OpenMPI/4.1.4-GCC-12.2.0
module load GROMACS/2021.5-foss-2021b

export GMX_MAXBACKUP=-1
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# 1) Energy minimization
gmx grompp -f ../../mdp/minim.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 1
mpirun gmx_mpi mdrun -v -deffnm em

# 2) NVT equilibration
gmx grompp -f ../../mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
mpirun gmx_mpi mdrun -deffnm nvt

# 3) NPT equilibration
gmx grompp -f ../../mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
mpirun gmx_mpi mdrun -deffnm npt

# 4) Production at 325 K
# Output name: mmc_{n_mmc}_{temp}
gmx grompp -f ../../mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -o mmc_{n_mmc}_{temp}.tpr -maxwarn 1
mpirun gmx_mpi mdrun -deffnm mmc_{n_mmc}_{temp}
"""

    for n in concentrations:
        folder_name = f"{n}_mmc"
        folder_path = os.path.join(conc_study_dir, folder_name)
        
        if not os.path.exists(folder_path):
            print(f"Skipping {folder_name} (directory not found)")
            continue
            
        script_content = template.format(n_mmc=n, temp=temperature)
        script_name = f"run_{n}_mmc.script"
        script_path = os.path.join(folder_path, script_name)
        
        with open(script_path, 'w', newline='\n') as f:
            f.write(script_content)
        
        print(f"Created {script_name} in {folder_name}")

if __name__ == "__main__":
    main()
