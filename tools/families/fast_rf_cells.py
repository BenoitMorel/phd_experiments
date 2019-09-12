import sys
import os
import math
import tempfile
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "print"))
sys.path.insert(0, os.path.join("tools", "families"))
import saved_metrics
import fam
import pickle
from aligned_printer import AlignedPrinter
import experiments as exp
import subprocess


def get_runs(datadir, run_tag):
  runs = []
  runs.append("true.true")
  if (run_tag == "all"):
    successful_runs = fam.get_successful_runs(datadir)
    for run in successful_runs:
      if (not "ultiple" in run and not "true" in run):
        runs.append(run)
  else:
    runs.append(run_tag)
  return runs
  
def get_trees(runs):
  trees = []
  for run in runs:
    trees.append(run + ".geneTree.newick")
  return trees

def write_rfdistance_input_files(families, trees, families_file, trees_file):
  with open(trees_file, "w") as writer:
    writer.write("\n".join(trees))
  with open(families_file, "w") as writer:
    for family in families:
      writer.write(family + " " + fam.get_gene_tree_dir(datadir, family) + "\n")

def run_rfdistance(families_file, trees_file, output_dir):
  command = []
  command.append(exp.rfdistance_exec)
  command.append(families_file)
  command.append(trees_file)
  command.append(output_dir)
  print(" ".join(command))
  subprocess.check_call(command)

def extract_rfdistance(rf_output_dir, families, runs):
  per_run_rf_vector = {}
  for run in runs: 
    per_run_rf_vector[run] = []
  for family in families:
    rf_file = os.path.join(rf_output_dir, family + ".rf")
    with open(rf_file) as reader:
      rf_distances = reader.readline()[:-2].split(" ")
      assert(len(rf_distances) == len(runs))
      for i in range(0, len(rf_distances)):
        per_run_rf_vector[runs[i]].append(float(rf_distances[i]))
  for run in runs:
    l = per_run_rf_vector[run]
    av = sum(l) / len(l)
    print(run + " " + str(av))

def analyze(datadir, run_tag, benched_run):
  temp_dir = tempfile.mkdtemp()#tempfile.TemporaryDirectory()
  analyze_dir_name = temp_dir#.name
  print("analyze directory " + analyze_dir_name)
  trees_file = os.path.join(analyze_dir_name, "trees.txt")
  families_file = os.path.join(analyze_dir_name, "families.txt")
  rf_output_dir = os.path.join(analyze_dir_name, "rfdistances")
  os.mkdir(rf_output_dir)

  families = fam.get_families_list(datadir)
  runs = get_runs(datadir, run_tag)
  print("Runs: " + str(runs))
  trees = get_trees(runs)
  write_rfdistance_input_files(families, trees, families_file, trees_file)
  run_rfdistance(families_file, trees_file, rf_output_dir)
  extract_rfdistance(rf_output_dir, families, runs)


if __name__ == '__main__':
  if (len(sys.argv) < 2):
    print("Syntax: families_dir run=all [benched_run]")
    exit(1)
  print(" ".join(sys.argv))
  datadir = sys.argv[1]
  run = sys.argv[2]
  benched_run = ""
  if (len(sys.argv) > 2):
    benched_run = sys.argv[2]
  analyze(datadir, run, benched_run)






