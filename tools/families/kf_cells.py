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
      if (not "ultiple" in run and not "true.true" in run):
        runs.append(run)
  else:
    runs.append(run_tag)
  return runs
  
def get_trees(runs):
  trees = []
  for run in runs:
    trees.append(run + ".geneTree.newick")
  return trees

def write_kfdistance_input_files(datadir, trees, trees_file):
  with open(trees_file, "w") as writer:
    writer.write("\n".join(trees))
    writer.write("\n")

def run_kfdistance(families_dir, trees_file, output_dir):
  command = []
  command.append("R")
  command.append("-f")
  command.append(exp.treedist_R_script)
  command.append("--args")
  command.append(families_dir)
  command.append(trees_file)
  command.append(output_dir)
  print(" ".join(command))
  subprocess.check_call(command)

def extract_kfdistance(kf_output_dir, families, runs):
  kf_cells = {}
  paired_runs = []
  for run in runs:
    paired_runs.append("true.true - " + run)
  for family in families:
    per_run_rkf = {}
    kf_cells[family] = per_run_rkf
    print("family " + str(family))
    kf_file = os.path.join(kf_output_dir, family)
    print(kf_file)
    with open(kf_file) as reader:
      kf_distances = reader.readline()[:-1].split(" ")
      print(kf_distances)
      assert(len(kf_distances) == len(runs))
      for i in range(0, len(kf_distances)):
        kf_cells[family][paired_runs[i]] = [float(kf_distances[i]), 1.0]
  return kf_cells

def save_kf_cells(datadir, kf_cells, rooted):
  output = fam.get_raw_kf_cells_file(datadir, rooted)
  pickle.dump(kf_cells, open(output, "wb"))

def load_kf_cells(datadir, rooted = False):
  return pickle.load(open(fam.get_raw_kf_cells_file(datadir, rooted), "rb"))

def get_run_key(m1, m2):
  return m1 + " - " + m2

def get_kf_to_true(cells, run_name):
  return cells[get_run_key(fam.get_run_name("true", "true"), run_name)]


def print_metrics(datadir, metric_dict, metric_name, benched_run):
  printer = AlignedPrinter()
  saved_metrics.save_dico(datadir, metric_dict, metric_name)
  for run_key in metric_dict:
    run = run_key.split(" - ")[1]
    suffix = ""
    if (benched_run == run):
      suffix += "\t <-- "
    printer.add("- " + run_key + ":",  str(metric_dict[run_key]) + suffix)
  printer.sort_right_float()
  printer.display()
  print("")


def export_metrics(datadir, benched_run, kf_cells, runs):
  total_rkf = {}
  families_number = len(kf_cells)
  run_keys = []
  for run in runs:
    run_keys.append("true.true - " + run)
  for run_key in run_keys:
    total_rkf[run_key] = 0.0
  for family in kf_cells:
    family_kf_cells = kf_cells[family]
    for key in family_kf_cells:
      total_rkf[key] += (family_kf_cells[key][0] / family_kf_cells[key][1]) 
  average_rkf = {}
  for key in run_keys:
    average_rkf[key] = total_rkf[key] / families_number
  print("Average KF:")
  print_metrics(datadir, average_rkf, "average_kf", benched_run)


def analyze(datadir, run_tag, benched_run = ""):
  temp_dir = tempfile.mkdtemp(dir = datadir)#tempfile.TemporaryDirectory()
  analyze_dir_name = temp_dir#.name
  print("analyze directory " + analyze_dir_name)
  trees_file = os.path.join(analyze_dir_name, "trees.txt")
  families_file = os.path.join(analyze_dir_name, "families.txt")
  kf_output_dir = os.path.join(analyze_dir_name, "kfdistances")
  os.mkdir(kf_output_dir)

  families_dir = fam.get_families_dir(datadir)
  families = fam.get_families_list(datadir)
  runs = get_runs(datadir, run_tag)
  print("Runs: " + str(runs))
  trees = get_trees(runs)
  write_kfdistance_input_files(datadir, trees, trees_file)
  
  run_kfdistance(families_dir, trees_file, kf_output_dir)
  kf_cells = extract_kfdistance(kf_output_dir, families, runs)
  if ("all" == run_tag):
    save_kf_cells(datadir, kf_cells, False)
  export_metrics(datadir, benched_run, kf_cells, runs) 

if __name__ == '__main__':
  if (len(sys.argv) < 3):
    print("Syntax: families_dir run=all [benched_run]")
    exit(1)
  print(" ".join(sys.argv))
  datadir = sys.argv[1]
  run = sys.argv[2]
  benched_run = ""
  if (len(sys.argv) > 3):
    benched_run = sys.argv[3]
  analyze(datadir, run, benched_run)






