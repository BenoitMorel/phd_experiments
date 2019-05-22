import sys
import os
import common
import pickle
sys.path.insert(0, os.path.join("tools", "families"))
from run_all import RunFilter

def run_single_experiment(dataset, do_generate, cores, run_filter = RunFilter()):
  if (do_generate):
    common.generate_dataset(dataset)
  common.run_reference_methods(dataset, cores, run_filter)
  #common.compute_likelihoods([dataset], cores)  
  common.run_all_analyzes([dataset])


if (__name__ == "__main__"): 
  if (len(sys.argv) < 4):
    print("Syntax: dataset do_generate cores [run_filter_pickled]")
    exit(1)

  dataset = sys.argv[1]
  do_generate = int(sys.argv[2])
  cores = int(sys.argv[3])
  run_filter = RunFilter()
  if (len(sys.argv) > 4):
    run_filter_file = sys.argv[4]
    print("run_filter pickle file: " + run_filter_file)
    run_filter = pickle.load(open(run_filter_file, "rb"))
    print("EXA_chains: " + str(run_filter.EXA_chains))

  run_single_experiment(dataset, do_generate, cores, run_filter)

