An MPI program that performs useless computations and MPI_Reduces. Used to experiment with some MPI schedulers.

Compilation:
mkdir build
cd build
cmake ..
make

Usage:
mpirun -np ranks fake_mpi_program number_of_trees number_of_sites

number of trees: number of iterations, each followed with a MPI_Reduce call.
number os sites: distributed among the cores
(there are no real trees/ sites in the program)



