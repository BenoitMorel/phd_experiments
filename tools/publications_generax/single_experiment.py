import sys
import os
import common
sys.path.insert(0, os.path.join("tools", "families"))
from run_all import RunFilter

def run_single_experiment(dataset, do_generate, cores):
  if (do_generate):
    common.generate_dataset(dataset)
  run_filter = RunFilter(raxml = False, pargenes = False, treerecs = False, phyldog = False)
  common.run_reference_methods(dataset, cores, run_filter)
  common.compute_likelihoods([dataset], cores)  
  common.run_all_analyzes([dataset])


if (__name__ == "__main__"): 
  if (len(sys.argv) != 4):
    print("Syntax: dataset do_generate cores")
    exit(1)

  dataset = sys.argv[1]
  do_generate = int(sys.argv[2])
  cores = int(sys.argv[3])
  run_single_experiment(dataset, do_generate, cores)
