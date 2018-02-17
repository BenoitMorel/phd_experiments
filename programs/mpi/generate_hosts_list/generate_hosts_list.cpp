#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <cstdlib>
#include <string>
#include <unordered_set>
using namespace std;

void master_main(int argc, char** argv)
{
  char *ranks = argv[1];
  string command =
    string("mpirun -np ")
    + string(ranks)
    + string(" /home/morelbt/github/phd_experiments/programs/mpi/generate_hosts_list/build/generate_hosts_list");
  system(command.c_str());
}

void printers_main(int argc, char** argv)
{
  MPI_Init(NULL, NULL);
  
  int name_len;
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(processor_name, &name_len);
  
  char *recieve_data = new char[MPI_MAX_PROCESSOR_NAME * world_size];

  MPI_Gather(processor_name,
      MPI_MAX_PROCESSOR_NAME,
      MPI_CHAR,
      recieve_data,
      MPI_MAX_PROCESSOR_NAME,
      MPI_CHAR,
      0,
      MPI_COMM_WORLD);
  if (!rank) {
    unordered_set<string> set;
    for (unsigned int i = 0; i < world_size; ++i) {
      set.insert(string(recieve_data + i * MPI_MAX_PROCESSOR_NAME));
    }
    for (auto str: set) {
      cout << str << endl;
    }
  }
  delete[] recieve_data;
  MPI_Finalize();
}

int main(int argc, char** argv) {
  if (argc == 1) {
    printers_main(argc, argv);
  } else if (argc == 2) {
    master_main(argc, argv);
  }
  return 0;
}

