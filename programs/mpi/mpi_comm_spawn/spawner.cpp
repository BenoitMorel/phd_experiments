/**
 * A fake program that performs fake computations
 * on fake trees and fake sites. It computes superfunctions
 * for each site on one given tree, and does a reduce on the
 * results. Then it iterates again on each tree. 
 * Each rank has a subset of sites.
 */

#include <mpi.h>
#include <iostream>
#include <chrono>
#include <vector>
#include <string>
#include <unistd.h>
using namespace std;


void spawn(const string &executable,
    const vector<string> &arguments,
    int spawns_number)
{
  char **argv = new char*[arguments.size() + 2];
  //argv[0] = (char*)"plop";
  for(unsigned int i = 0; i < arguments.size(); ++i)
    argv[i] = (char*)arguments[i].c_str();
  argv[arguments.size()] = 0;

  MPI_Comm intercomm;
  MPI_Comm_spawn((char*)executable.c_str(), argv, spawns_number,  
          MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm,  
          MPI_ERRCODES_IGNORE);
  delete[] argv;

}

void spawner()
{
  cout << "I am the spawner" << endl;
  // string executable = "../../fake_mpi_program/build/fake_mpi_program";
  string executable = "../script_wrap.sh";
  vector<string> arguments;
  arguments.push_back("1000");
  arguments.push_back("100000");
  
  int spawns_number = 3;
  MPI_Comm intercomm;
  for (unsigned int iter = 0; iter < 3; ++iter) { 
    std::cout << "spawner: iter " << iter << std::endl;
    for (unsigned int i = 0; i < 2; ++i) {
      spawn(executable, arguments, spawns_number);
    }
    usleep(10 * 1000 * 1000);
  }
  cout << "End of spawner" << std::endl;
}

void spawned()
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  cout << "I am the spawned " << rank << endl;
  usleep(3000000);
  cout << "end" << endl;
}


int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);        
  MPI_Comm parent;
  MPI_Comm_get_parent(&parent); 
  if (parent == MPI_COMM_NULL) {
    spawner();
  } else {
    spawned();
}
  /*
  auto start = chrono::system_clock::now();
 
  if (argc != 3) {
    cerr << "Invalid syntax." << endl;
    cerr << "fake_mpi_program <trees_number> <sites_number>" << endl;
    exit(0);
  }

  int trees = atoi(argv[1]);
  int sites = atoi(argv[2]);
  fake_mpi_program(trees, sites);

  auto end = chrono::system_clock::now();
  auto elapsed = chrono::duration_cast<chrono::milliseconds>
                                 (end-start).count();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    cout << "Elapsed time " << elapsed << "ms" << endl; 
  
  */  
  MPI_Finalize();              
  return 0;
}


