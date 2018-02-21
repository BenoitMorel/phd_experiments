#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>

int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get the name of the processor
  int name_len;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  char hostname[HOST_NAME_MAX];
  MPI_Get_processor_name(processor_name, &name_len);
  gethostname(hostname, HOST_NAME_MAX);
  

  // Print off a hello world message
  printf("%s\n", processor_name);
  printf("%s\n", hostname);

  // Finalize the MPI environment.
  MPI_Finalize();
}

