#!/bin/bash
echo "Starting sequential runs in background. Check run_all.log"
nohup bash -c '
set -e
echo "Processing 2_mmc..."
cd 2_mmc
gmx grompp -f ../mdp/minim.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 1 && gmx mdrun -v -deffnm em -nt 32\ngmx grompp -f ../mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1 && gmx mdrun -deffnm nvt -nt 32\ngmx grompp -f ../mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1 && gmx mdrun -deffnm npt -nt 32\ngmx grompp -f ../mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -o mmc_2_325K.tpr -maxwarn 1 && gmx mdrun -deffnm mmc_2_325K -nt 32\ncd ..\necho "Processing 4_mmc..."
cd 4_mmc
gmx grompp -f ../mdp/minim.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 1 && gmx mdrun -v -deffnm em -nt 32\ngmx grompp -f ../mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1 && gmx mdrun -deffnm nvt -nt 32\ngmx grompp -f ../mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1 && gmx mdrun -deffnm npt -nt 32\ngmx grompp -f ../mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -o mmc_4_325K.tpr -maxwarn 1 && gmx mdrun -deffnm mmc_4_325K -nt 32\ncd ..\necho "Processing 8_mmc..."
cd 8_mmc
gmx grompp -f ../mdp/minim.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 1 && gmx mdrun -v -deffnm em -nt 32\ngmx grompp -f ../mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1 && gmx mdrun -deffnm nvt -nt 32\ngmx grompp -f ../mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1 && gmx mdrun -deffnm npt -nt 32\ngmx grompp -f ../mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -o mmc_8_325K.tpr -maxwarn 1 && gmx mdrun -deffnm mmc_8_325K -nt 32\ncd ..\necho "Processing 12_mmc..."
cd 12_mmc
gmx grompp -f ../mdp/minim.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 1 && gmx mdrun -v -deffnm em -nt 32\ngmx grompp -f ../mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1 && gmx mdrun -deffnm nvt -nt 32\ngmx grompp -f ../mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1 && gmx mdrun -deffnm npt -nt 32\ngmx grompp -f ../mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -o mmc_12_325K.tpr -maxwarn 1 && gmx mdrun -deffnm mmc_12_325K -nt 32\ncd ..\necho "Processing 16_mmc..."
cd 16_mmc
gmx grompp -f ../mdp/minim.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 1 && gmx mdrun -v -deffnm em -nt 32\ngmx grompp -f ../mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1 && gmx mdrun -deffnm nvt -nt 32\ngmx grompp -f ../mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1 && gmx mdrun -deffnm npt -nt 32\ngmx grompp -f ../mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -o mmc_16_325K.tpr -maxwarn 1 && gmx mdrun -deffnm mmc_16_325K -nt 32\ncd ..\n' > run_all.log 2>&1 &
