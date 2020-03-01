This folder contains some of the files that I used for the NAMD tutorial

https://www.youtube.com/watch?v=yKAMRMigmvU

I have included the
1) Topology file
2) Parameter file
3) PDB file

Instructions are given in the video to create the NAMD script used for energy minimization, equilibration and NVT simulations.

Example 1:
1. Load 5eui.pdb into VMD
2. Create PSF file using 'top_all27_prot_lipid_na.inp'
3. Step 2 should give you similar files like original.pdb and original.psf
4. use the par_all27_prot_lipid_na.inp in the conf file (runme.namd) to run simulation.
5. Note that the parameters in runme.namd are default and may not correct, especially for PBC box dimensions.

Hope this helps!
