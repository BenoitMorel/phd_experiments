# Common variables used in the scripts

import os
import datetime
import sys

root = "./"

# scripts
scripts_root = os.path.join(root, "scripts")

# datasets
datasets_root = os.path.join(root, "datasets")

# results
results_root = os.path.join(root, "results")

# externals
treerecs_root = os.path.join("..", "Treerecs")

# utils
def write_results_info(resultsdir, msg):
  filename = os.path.join(resultsdir, "results.info")
  with open(filename, "w") as writer:
      writer.write("Experiment started on "+ str(datetime.datetime.now()) + "\n")
      writer.write("Command: " + ' '.join(sys.argv + "\n"))
      writer.write(msg)
  



