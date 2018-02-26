#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#include <chrono>
#include <limits.h>
#include <cstdlib>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/wait.h>
using namespace std;


string myself = "/home/morelbt/github/phd_experiments/programs/mpi/mpi_and_exec/build/mpi_and_exec";
string expansiveExec = "/home/morelbt/github/phd_experiments/programs/mpi/fake_mpi_program/build/fake_mpi_program 500 50000";

void remove_domain(string &str)
{
  str = str.substr(0, str.find(".", 0));
}

void master_main(int argc, char** argv)
{

  const char* executable = myself.c_str();
  const char * arguments[3];
  arguments[0] = executable;
  arguments[1] = "slave";
  arguments[2] = (const char*)0;
  cout << "This experiment is a failure because:" << endl;
  cout << "- execv does not return in case of success" << endl;
  cout << "- execv does not allow to pass a MPI context to the new process" << endl;
  cout << "exe 1 " << execv(executable, (char* const*)arguments) << endl;

  cout << "hey execv finished" << endl;
}

void slave_main(int argc, char** argv)
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
    map<string, int> nodes;;
    for (unsigned int i = 0; i < world_size; ++i) {
      string node(recieve_data + i * MPI_MAX_PROCESSOR_NAME);
      remove_domain(node);
      int size = 0;
      if (nodes.find(node) != nodes.end()) {
        size = nodes.at(node);
      }
      nodes[node] = ++size;
    }
    for (auto pair: nodes) {
      cout << pair.first << " slots = " << pair.second << endl;
    }
  }
  delete[] recieve_data;
  MPI_Finalize();
}

int main(int argc, char** argv) {
  if (argc == 1) {
    cerr << "invalid syntax" << endl;
    return 1;
  }
  if (string(argv[1]) == string("slave")) {
    slave_main(argc, argv);
  } else if (string(argv[1]) == string("master")) {
    master_main(argc, argv);
  }
  return 0;
}

