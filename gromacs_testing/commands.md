<!--
Martini3 Fibrinogen Simulation Workflow
Author: Sumith Yesudasan
Purpose: Fully reproducible documentation of CG simulation workflow
-->

# Martini3 Fibrinogen Simulation Workflow (GROMACS + Ubuntu)

This document describes a **complete workflow for running Martini3 coarse-grained simulations of fibrinogen using GROMACS**.

The steps include:

1. Environment setup  
2. Protein coarse-graining  
3. System building  
4. Energy minimization  
5. Equilibration  
6. Production simulation  
7. Trajectory processing  

---

# 1. System Requirements

Recommended system:

- Ubuntu 20.04 / 22.04  
- Python 3.10  
- Conda / Miniconda  
- GROMACS 2023+  
- Martini3 force field  
- martinize2  
- insane  

Optional visualization software:

- VMD  
- OVITO  
- PyMOL  

---

# 2. Create Simulation Environment

Create a dedicated Conda environment.

```bash
conda create -n martini3 python=3.10
conda activate martini3
```

Install simulation tools:

```bash
conda install -c conda-forge gromacs
pip install martinize2 insane
```

## Verify Installation

```bash
gmx --version
insane -h
martinize2 -h
```

### Expected Output

```
GROMACS version information
INSANE command help
martinize2 command help
```

---

# 3. Prepare Atomistic Protein Structure

Input file:

```
fibrinogen.pdb
```

Check structure:

```bash
head fibrinogen.pdb
```

### Expected Output

The file should contain lines beginning with:

```
ATOM
HETATM
```

---

# 4. Convert Atomistic Protein to Martini3 CG

Run martinize2:

```bash
martinize2 -f fibrinogen.pdb -x fibrinogen_cg.pdb -o molecule.itp -ff martini3001 -water W
```

### Expected Output Files

```
fibrinogen_cg.pdb
molecule_0.itp
molecule_1.itp
molecule_2.itp
molecule_3.itp
molecule_4.itp
```

These files represent:

| File | Description |
|-----|-------------|
| fibrinogen_cg.pdb | CG coordinates |
| molecule_*.itp | topology for protein chains |

---

# 5. Build Simulation Box and Solvate System

Use **INSANE** to create the solvated simulation system.

```bash
insane -f fibrinogen_cg.pdb -o CG.gro -p system.top -sol W -salt 0.15 -d 8 -pbc cubic
```

### Explanation

| Parameter | Meaning |
|-----------|--------|
| -f | CG protein structure |
| -o | output coordinate file |
| -p | topology file |
| -sol W | Martini water beads |
| -salt 0.15 | physiological salt |
| -d 8 | protein-box distance |
| -pbc cubic | cubic periodic box |

### Expected Output Files

```
CG.gro
system.top
```

---

# 6. Check Topology File

Open topology:

```bash
nano system.top
```

Expected includes:

```
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"
```

Expected molecule definitions:

```
[molecules]
Protein_chain_A    1
Protein_chain_B    1
Protein_chain_C    1
W                  XXXX
NA                 XX
CL                 XX
```

---

# 7. Energy Minimization

Prepare simulation input:

```bash
gmx grompp -f minimization.mdp -c CG.gro -r CG.gro -p system.top -o minimization.tpr
```

Run minimization:

```bash
gmx mdrun -deffnm minimization -v
```

### Expected Output Files

```
minimization.tpr
minimization.gro
minimization.edr
minimization.log
```

---

# 8. Equilibration Stage 1

```bash
gmx grompp -f equil1.mdp -c minimization.gro -r minimization.gro -p system.top -o equil1.tpr -DPOSRES
```

Run:

```bash
gmx mdrun -deffnm equil1
```

Expected output:

```
equil1.tpr
equil1.gro
equil1.edr
equil1.log
equil1.xtc
```

---

# 9. Equilibration Stage 2

```bash
gmx grompp -f equil2.mdp -c equil1.gro -r equil1.gro -p system.top -o equil2.tpr -DPOSRES
```

Run:

```bash
gmx mdrun -deffnm equil2
```

---

# 10. Equilibration Stage 3

```bash
gmx grompp -f equil3.mdp -c equil2.gro -p system.top -o equil3.tpr
```

Run:

```bash
gmx mdrun -deffnm equil3
```

---

# 11. Production Simulation

```bash
gmx grompp -f production.mdp -c equil3.gro -p system.top -o production.tpr
```

Run:

```bash
gmx mdrun -deffnm production
```

Expected output:

```
production.tpr
production.xtc
production.edr
production.log
production.gro
```

---

# 12. Trajectory Processing

```bash
gmx trjconv -f production.xtc -s production.tpr -o centered.xtc -center -pbc mol
```

---

# 13. Visualization

```bash
vmd centered.xtc
```

or

```bash
ovito centered.xtc
```

---

# 14. Useful Ubuntu Commands

```bash
du -sh *
htop
ls -lh
rm step*.pdb
```

---

# 15. Recommended Project Directory Structure

```
fibrinogen-martini/
│
├── README.md
├── input/
│   └── fibrinogen.pdb
│
├── mdp/
│   ├── minimization.mdp
│   ├── equil1.mdp
│   ├── equil2.mdp
│   ├── equil3.mdp
│   └── production.mdp
│
├── topology/
│   ├── molecule_0.itp
│   ├── molecule_1.itp
│   ├── molecule_2.itp
│
├── run/
│   ├── CG.gro
│   ├── system.top
│   └── trajectories/
│
└── analysis/
```

---

# 16. Troubleshooting

## LINCS Errors

Possible causes:

- timestep too large  
- poor minimization  
- solvent overlaps  

Solutions:

- reduce timestep  
- run longer minimization  
- apply position restraints  
