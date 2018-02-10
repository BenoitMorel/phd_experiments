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

using namespace std;

double superfunction(unsigned int s, const vector<double> &v)
{
  double res = 0;
  for (auto d: v)
    res += d * s + s * s;
  return res;
}

void fake_mpi_program(int trees, int sites) 
{

  int rank;
  int size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  if (rank == 0) {
    cout << "Fake computation on " << trees << " trees " <<
      "and " << sites << " sites with " << size << " ranks."<< endl;

  }
  
  vector<double> bidon_input;
  for (unsigned int i = 0; i < 100; ++i) {
    bidon_input.push_back(i * 5 + i);
  }
  double globalTreell = 0.0;
  int start = (sites / size) * rank;
  int end = min(sites, start + sites / size);
  for (unsigned int tree = 0; tree < trees; ++tree) {
    double localTreell = 0.0;
    for (unsigned int s = start; s < end; ++s) {
      localTreell += superfunction(s, bidon_input);      
    }
    double temp;
    MPI_Reduce(&localTreell, &temp, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
    globalTreell += temp;
  }

  if (rank == 0) {
    cout << "result: " <<  globalTreell << endl;
  }

}

int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);        
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
  cout << "Elapsed time " << elapsed << "ms" << endl; 
  MPI_Finalize();              
  return 0;
}


