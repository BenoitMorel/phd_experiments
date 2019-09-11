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

def analyze(datadir, method, benched_run):
  temp_dir = tempfile.mkdtemp()#tempfile.TemporaryDirectory()
  analyze_dir_name = temp_dir#.name
  print("analyze directory " + analyze_dir_name)
  trees_file = os.path.join(analyze_dir_name, "trees.txt")
  families_file = os.path.join(analyze_dir_name, "families.txt")
  output_dir = os.path.join(analyze_dir_name, "rfdistances")
  os.mkdir(output_dir)

  families = fam.get_families_list(datadir)
  methods = []
  if (method == "all"):
    methods = fam.get_successful_runs(datadir)
  else:
    methods.append(method)
  trees = []
  trees.append("true.true.geneTree.newick")
  for method in methods:
    if (method != "true.true"):
      trees.append(method + ".geneTree.newick")
  with open(trees_file, "w") as writer:
    writer.write("\n".join(trees))
  with open(families_file, "w") as writer:
    for family in families:
      writer.write(fam.get_gene_tree_dir(datadir, family) + "\n")
  command = []
  command.append(exp.rfdistance_exec)
  command.append(families_file)
  command.append(trees_file)
  command.append(output_dir)
  subprocess.check_call(command)

if __name__ == '__main__':
  if (len(sys.argv) < 2):
    print("Syntax: families_dir method=all [benched_run]")
    exit(1)
  print(" ".join(sys.argv))
  datadir = sys.argv[1]
  method = sys.argv[2]
  benched_run = ""
  if (len(sys.argv) > 2):
    benched_run = sys.argv[2]
  analyze(datadir, method, benched_run)






