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


string myself = "/home/morelbt/github/phd_experiments/programs/mpi/generate_hosts_list/build/generate_hosts_list";
string printer_output = "/home/morelbt/github/phd_experiments/programs/mpi/generate_hosts_list/printer_output.txt";
string expansiveExec = "/home/morelbt/github/phd_experiments/programs/mpi/fake_mpi_program/build/fake_mpi_program 500 50000";

void remove_domain(string &str)
{
  str = str.substr(0, str.find(".", 0));
}

int forkAndGetPid(const string & command)
{

  int pid = fork();
  if (!pid) {
    system(command.c_str());
    exit(0);
  } else {
    return pid;
  }
}

void master_main(int argc, char** argv)
{
  using milli = std::chrono::milliseconds;
  auto start = std::chrono::high_resolution_clock::now();
  char *ranks = argv[1];
  string commandPrinters =
    string("mpirun -np ")
    + string(ranks)
    + string(" ")
    + myself
    + string (" printer ")
    + string("> ")
    + printer_output;
  system(commandPrinters.c_str());
  ifstream hostfileReader(printer_output);
  vector<pair<string, int> > nodes;
  while (!hostfileReader.eof()) {
    string host;
    string temp;
    int slots;
    hostfileReader >> host;
    if (hostfileReader.eof())
      break;
    hostfileReader >> temp;
    hostfileReader >> temp;
    hostfileReader >> slots;
    nodes.push_back(pair<string, int>(host, slots));
  }
  vector<int> pids;
  for (auto pair: nodes) { 
    cout << pair.first << " " << pair.second << endl;
    string commandRun = 
      string ("mpirun ")
      
      + string("--hostfile ")
      + printer_output
      + string(" -host ")
      + pair.first 
      //+ string(" -display-allocation ")
      + string(" -np ")
      + to_string(pair.second)
      + string(" ")
      + expansiveExec;
    cout << "execute " << commandRun << endl;
    pids.push_back(forkAndGetPid(commandRun.c_str()));
  }
  auto finish = std::chrono::high_resolution_clock::now();
  std::cout << "before wait "
    << std::chrono::duration_cast<milli>(finish - start).count()
    << " milliseconds\n";
  for (auto pid: pids)
    waitpid(pid, 0, 0);
  //system("wait");
  finish = std::chrono::high_resolution_clock::now();
  std::cout << "total "
    << std::chrono::duration_cast<milli>(finish - start).count()
    << " milliseconds\n";
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
  if (string(argv[1]) == string("printer")) {
    printers_main(argc, argv);
  } else {
    master_main(argc, argv);
  }
  return 0;
}

