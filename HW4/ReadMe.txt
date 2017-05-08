In this assignment,I improved the performance of two applications, Conjugate Gradient and Multi-Grid .we know that Conjugate Gradient and Multi-Grid are two kernel applications in NAS Parallel  Benchmarks(NPB).  NPB is a small set of programs designed to help evaluate the performance of parallel supercomputers.Downloaded the benchmark ,that can be foundhttp://www.nas.nasa.gov/publications/npb.html.

              
I improved the performance of CG and MG by making changes in OpenMp pragmas,added few OpenMP constructs and removed few OpenMP constructs ,and added environment varibles in the code.

I have included make file in Parallel code of both CG and MG.

The modified version of code can be build and run CG AND MG (Class)by  folwing instructions:


1.Conjugate Gradient(CG)

$ tar xfz HW4.tar.gz

$ cd HW4/ParallelCode/NPB3.3-OMP-C/

$ make cg CLASS=C

$./bin/cg.C.x




2.Multi-Grid(MG)



$ tar xfz HW4.tar.gz

$ cd HW4/ParallelCode/NPB3.3-OMP-C/

$ make mg CLASS=C

$ ./bin/mg.C.x
