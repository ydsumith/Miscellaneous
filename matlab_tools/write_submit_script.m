function write_submit_script(runfile)
outfile_temp = fullfile(runfile,'submit.sh');

fid = fopen(outfile_temp, 'w');
fprintf(fid,'#PBS -S /bin/bash\n');
fprintf(fid,'#PBS -j oe\n');
fprintf(fid,'#PBS -q batch\n');
fprintf(fid,'#PBS -N %s\n', runfile);
fprintf(fid,'#PBS -l nodes=1:ppn=6:rdanode\n');
fprintf(fid,'#PBS -l walltime=48:00:00\n');
fprintf(fid,'#PBS -l vmem=25g\n');
fprintf(fid,'cd $PBS_O_WORKDIR\n');
fprintf(fid,'module load openmpi/1.8.3/gcc/4.4.7\n');
fprintf(fid,'mpirun -np $PBS_NP /lustre1/sy49293/lammps/lammps-09MAR2018/src/lmp_mpi  < conf.in > %s.out\n',runfile);
fclose(fid);

end