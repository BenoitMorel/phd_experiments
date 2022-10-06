import os
import sys
sys.path.insert(0, 'tools/plotters')
import plot_line
import get_raxmlruntimes
 

def get_treecombine_runtime(rundir):
  logs = os.path.join(rundir, "logs.txt")
  print(logs)
  lb = 0.0
  total = 0.0
  cores = 0
  line_number = 0
  for line in open(logs).readlines():
    if (line_number == 3):
      assert("mpi-scheduler" in line)
      cores = float(line.split()[2])
    if (line.startswith("Finished running commands.")):
      total = float(line.split()[-1][:-1])
    if (line.startswith("Load balance ratio:")):
      lb = float(line.split()[-1])
    line_number += 1
  print(lb)
  print(total)
  print(cores)
  print(total * cores * lb)
  return total * cores * lb

def plot(datadir, totaltreenumber, treenumbers, output):
  
  yraxml = []
  ytreecombine = []
   
  treebase = "/hits/fast/cme/schmidja/difficulty_training_data/TreeBase_results/AA/"
  print("Gathering raxml runtimes...")
  raxml_runtime = get_raxmlruntimes.get_raxmlruntimes(datadir, treebase)
  print("Gathering treecombine runtimes...")
  for n in treenumbers:
    yraxml.append(float(raxml_runtime) / (float(totaltreenumber) / float(n)))
    treecombine_run_dir = os.path.join(datadir, "runs", "LG+G", "treecombination" + str(totaltreenumber) + "_s" + str(n))
    treecombine_runtime = get_treecombine_runtime(treecombine_run_dir)
    ytreecombine.append(treecombine_runtime)
  yvalues = [yraxml, ytreecombine]
  line_captions = ["RAxML-NG", "TreeCombine"]
  plot_line.plot_line(treenumbers, yvalues, "TreeCombine VS RAxML-NG",  "Starting trees", "Runtime (s)", output, line_captions)


if (__name__ == "__main__"):
  if (len(sys.argv) < 4):
    print("Syntax python " + os.path.basename(__file__) + " datadir totaltreenumber output [treenumber1, treenumber2, ...]")
    sys.exit(1)
  datadir = sys.argv[1]
  totaltreenumber = int(sys.argv[2])
  output = sys.argv[3]
  treenumbers = []
  for n in sys.argv[4:]:
    treenumbers.append(int(n))

  plot(datadir, totaltreenumber, treenumbers, output)

