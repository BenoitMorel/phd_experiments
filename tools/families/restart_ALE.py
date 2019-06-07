import os
import sys
sys.path.insert(0, 'scripts')
import experiments as exp

import run_ALE

def launch(datadir, cluster, cores):
  dataset = os.path.basename(os.path.normpath(datadir)) 
  command = ["python"]
  command.extend(sys.argv)
  command.append("--exprun")
  resultsdir = os.path.join("RestartAle", dataset)
  resultsdir = exp.create_result_dir(resultsdir, [])
  submit_path = os.path.join(resultsdir, "submit.sh")
  exp.submit(submit_path, " ".join(command), cores, cluster) 


if (__name__ == "__main__"):
  is_run = ("--exprun" in sys.argv)
  if (len(sys.argv) < 4):
    print("Syntax: datadir cores cluster")
    exit(1)
  datadir = sys.argv[1]
  cores = int(sys.argv[2])
  if (is_run):
    run_ALE.restart_exabayes_and_ALE(datadir, cores)
  else:
    cluster = sys.argv[3]
    launch(datadir, cluster, cores)
