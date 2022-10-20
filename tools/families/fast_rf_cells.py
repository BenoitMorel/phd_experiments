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
      if (not "poly" in run and not "ultiple" in run and not "true.true" in run and not "disco" in run):
        runs.append(run)
  else:
    runs.append(run_tag)
  return runs
  
def get_trees(runs):
  trees = []
  for run in runs:
    trees.append(run + ".geneTree.newick")
  return trees

def write_rfdistance_input_files(datadir, families, trees, families_file, trees_file):
  with open(trees_file, "w") as writer:
    writer.write("\n".join(trees))
  with open(families_file, "w") as writer:
    for family in families:
      writer.write(family + " " + fam.get_gene_tree_dir(datadir, family) + "\n")

def run_rfdistance(families_file, trees_file, output_dir, cores):
  command = []
  command.append("mpiexec")
  command.append("-np")
  command.append(str(cores))
  command.append(exp.rfdistance_exec)
  command.append(families_file)
  command.append(trees_file)
  command.append(output_dir)
  print(" ".join(command))
  subprocess.check_call(command)

def extract_rfdistance(rf_output_dir, families, runs):
  rf_cells = {}
  paired_runs = []
  for run in runs:
    paired_runs.append("true.true - " + run)
  for family in families:
    per_run_rrf = {}
    rf_cells[family] = per_run_rrf
    rf_file = os.path.join(rf_output_dir, family + ".rf")
    with open(rf_file) as reader:
      max_rf_distances = float(reader.readline()[:-1])
      rf_distances = reader.readline()[:-1].split(" ")
      assert(len(rf_distances) == len(runs))
      for i in range(0, len(rf_distances)):
        rf_cells[family][paired_runs[i]] = [float(rf_distances[i]), float(max_rf_distances)]
  return rf_cells

def save_rf_cells(datadir, rf_cells, rooted):
  output = fam.get_raw_rf_cells_file(datadir, rooted)
  pickle.dump(rf_cells, open(output, "wb"))

def load_rf_cells(datadir, rooted = False):
  return pickle.load(open(fam.get_raw_rf_cells_file(datadir, rooted), "rb"))

def print_metrics(datadir, metric_dict, metric_name, benched_run):
  printer = AlignedPrinter()
  saved_metrics.save_dico(datadir, metric_dict, metric_name)
  for run_key in metric_dict:
    run = run_key #.split(" - ")[1]
    suffix = ""
    if (benched_run == run):
      suffix += "\t <-- "
    printer.add("- " + run_key + ":",  str(metric_dict[run_key]) + suffix)
  printer.sort_right_float()
  printer.display()
  print("")


def export_metrics(datadir, benched_run, rf_cells, runs):
  total_rrf = {}
  valid_trees = {}
  families_number = len(rf_cells)
  run_keys = []
  for run in runs:
    run_keys.append("true.true - " + run)
  for run_key in run_keys:
    total_rrf[run_key] = 0.0
    valid_trees[run_key] = 0.0
  for family in rf_cells:
    family_rf_cells = rf_cells[family]
    for key in family_rf_cells:
      is_valid = family_rf_cells[key][0] >= 0.0
      if (is_valid):
        total_rrf[key] += (family_rf_cells[key][0] / family_rf_cells[key][1])
        valid_trees[key] += 1.0
  average_rrf = {}
  
  for key in run_keys:
    invalid = families_number - valid_trees[key]
    if (invalid > 0):
      print("Warning: found " + str(invalid) + " invalid trees for method " + key)
    average_rrf[key.split(" - ")[1]] = total_rrf[key] / valid_trees[key]
  
  print("Average relative RF:")
  print_metrics(datadir, average_rrf, "average_rrf", benched_run)


def analyze(datadir, run_tag, cores, benched_run = ""):
  temp_dir = tempfile.mkdtemp(dir = datadir)#tempfile.TemporaryDirectory()
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
  write_rfdistance_input_files(datadir, families, trees, families_file, trees_file)
  run_rfdistance(families_file, trees_file, rf_output_dir, cores)
  rf_cells = extract_rfdistance(rf_output_dir, families, runs)
  if ("all" == run_tag):
    save_rf_cells(datadir, rf_cells, False)
  export_metrics(datadir, benched_run, rf_cells, runs) 

if __name__ == '__main__':
  if (len(sys.argv) < 4):
    print("Syntax: families_dir run=all cores [benched_run]")
    exit(1)
  print(" ".join(sys.argv))
  datadir = sys.argv[1]
  run = sys.argv[2]
  cores = int(sys.argv[3])
  benched_run = ""
  if (len(sys.argv) > 4):
    benched_run = sys.argv[4]
  analyze(datadir, run, cores, benched_run)






