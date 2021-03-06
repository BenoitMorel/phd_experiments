import time
import sys
import os
import fam_data
import pickle
sys.path.insert(0, os.path.join("tools", "families"))
from run_all import RunFilter

def run_single_experiment(dataset, subst_model, do_generate, cores, run_filter = RunFilter()):
  start = time.time()
  try:
    if (do_generate):
      fam_data.generate_dataset(dataset)
    run_filter.run_reference_methods(dataset, subst_model, cores)
  finally:
    elapsed = time.time() - start
    print("End of single experiment. Elapsed time: " + str(elapsed) + "s")

if (__name__ == "__main__"): 
  if (len(sys.argv) < 5):
    print("Syntax: dataset do_generate cores [run_filter_pickled]")
    exit(1)

  print("start single experiment " + " ".join(sys.argv))
  dataset = sys.argv[1]
  subst_model = sys.argv[2]
  do_generate = int(sys.argv[3])
  cores = int(sys.argv[4])
  run_filter = RunFilter()
  if (len(sys.argv) > 5):
    run_filter_file = sys.argv[5]
    print("run_filter pickle file: " + run_filter_file)
    run_filter = pickle.load(open(run_filter_file, "rb"))

  run_filter.run_single_experiment(dataset, subst_model, do_generate, cores)

