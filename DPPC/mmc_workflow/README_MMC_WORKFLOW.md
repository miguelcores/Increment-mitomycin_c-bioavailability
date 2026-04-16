MMC Insertion Workflow for 256-DPPC Bilayer (325 K, 1 bar)

Files generated in this folder:
- mmc.itp: united-atom MMC topology (renamed moleculetype MMC)
- topol_mmc.top: topology template with MMC include and molecule placeholders
- prepare_mmc_system.py: water replacement + MMC insertion with overlap rejection
- validate_mmc_system.py: checks for molecule counts, net charge, and estimated density
- concentration_to_counts.py: converts target mM values to DPPC/MMC/SOL molecule counts in this box
- mdp/*.mdp: EM/NVT/NPT/production settings at 325 K and 1 bar
- run_mmc_workflow.sh: end-to-end command sequence

How to build one concentration:
1) Convert concentration target(s) to molecule counts for this box
   python3 concentration_to_counts.py --mM 5 10 20 50

2) Build the inserted coordinates and auto-updated topology
   python3 prepare_mmc_system.py --n-mmc 32 --seed 325 \
      --output-gro system_32mmc.gro --output-top system_32mmc.top

3) Validate counts/charge/density estimate
   python3 validate_mmc_system.py --gro system_32mmc.gro --top system_32mmc.top

4) Run GROMACS stages (EM -> NVT -> NPT -> MD)
   ./run_mmc_workflow.sh 32 325

Molecule counting logic:
- DPPC is fixed at 256 molecules.
- Each inserted MMC removes waters-per-mmc solvent residues (default 2).
- SOL count in the output topology is automatically recomputed.

Overlap control:
- MMC is randomly rotated and translated to selected water oxygen anchors.
- Trial poses are accepted only when all inserted atoms are at least min-dist from existing atoms under periodic boundary conditions.
- If insertion fails, increase --waters-per-mmc or reduce --min-dist modestly.

Notes on force-field compatibility:
- The MMC file is ATB united-atom based on GROMOS 54A7 and has net charge 0.
- Non-bonded mixing remains controlled by forcefield.itp (Lorentz-Berthelot/comb-rule from the selected GROMOS lipid force field).
