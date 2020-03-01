Use the conf.gro and compile it using 
>>g++ *.cpp -fopenmp -o test.out

gromacs files can be run using 
>>grompp -v
>>mdrun -v

for charges
>>grompp -n index.ndx -v
>>mdrun -v

The results are uploaded in txt files.
It doesnt match the force values with gromacs
