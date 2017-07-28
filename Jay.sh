mpirun -npernode 1 ./mpifinal 2000 output.xt

qsub -cwd -pe mpich 4 jay.sh

