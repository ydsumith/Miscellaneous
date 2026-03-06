gmx grompp -f minimization.mdp -c CG.gro -r CG.gro -p system.top -o minimization.tpr
gmx mdrun -deffnm minimization -v
