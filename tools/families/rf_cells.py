import sys
import os
from ete3 import Tree
import numpy
import math
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "families"))
from read_tree import read_tree
from read_tree import read_trees_list
from rf_distance import ete3_rf
from rf_distance import get_relative_rf
from rf_distance import get_rf
import rf_distance
import concurrent.futures
import saved_metrics
import fam
import pickle

class AlignedPrinter:
  def __init__(self):
    self._lefts = []
    self._rights = []

  def add(self, left, right):
    self._lefts.append(left)
    self._rights.append(right)

  def sort_right_float(self):
    floating_array = []
    for elem in self._rights:
      floating_array.append(elem.split()[0])
    self._lefts = [x for _,x in sorted(zip(floating_array, self._lefts))]
    self._rights = [x for _,x in sorted(zip(floating_array, self._rights))]

  def display(self):
    max_chars = 0
    for l in self._lefts:
      max_chars = max(max_chars, len(l))
    for i in range(0, len(self._lefts)):
      l = self._lefts[i]
      r = self._rights[i]
      to_print = l + " " * (max_chars - len(l) + 1) + r
      print(to_print)

def get_gene_trees_list(datadir, family, run):
  tree_path = fam.build_gene_tree_path_from_run(datadir, family, run)
  return read_trees_list(tree_path)

def read_trees_for_family(datadir, family, runs):
  family_trees = {}
  for run in runs:
      family_trees[run] = get_gene_trees_list(datadir, family, run)
  return family_trees

def read_all_trees(datadir, runs):
  trees = {}
  families = fam.get_families_list(datadir)
  for family in families:
    trees[family] = read_trees_for_family(datadir, family, runs)
  return trees

def get_invalid_runs(trees, runs):
  invalid_runs = []
  for run in runs:
    for family in trees:
      if (not run in trees[family]):
        print("invalid run: " + run)
        invalid_runs.append(run)
        break
  return invalid_runs

def get_runs_to_compare(runs):
  runs_to_compare = []
  for run in runs:
    runs_to_compare.append(("true.true", run))
  return runs_to_compare

def get_run_key(m1, m2):
  return m1 + " - " + m2

def compute_all_rf_cells(trees, runs_to_compare):
  rf_cells = {}
  for family in trees:
    family_trees = trees[family]
    rf_cells[family] = {}
    for run1, run2 in runs_to_compare:
      key = get_run_key(run1, run2)
      trees1 = family_trees[run1]
      trees2 = family_trees[run2]
      rf_cells[family][key] = rf_distance.ete3_average_rf_from_list(trees1, trees2) 
  return rf_cells

def save_rf_cells(datadir, rf_cells):
  output = fam.get_raw_rf_cells_file(datadir)
  pickle.dump(rf_cells, open(output, "wb"))

def load_rf_cells(datadir):
  return pickle.load(open(fam.get_raw_rf_cells_file(datadir), "rb"))

def get_rf_to_true(cells, run_name):
  return cells[get_run_key(fam.get_run_name("true", "true"), run_name)]

def compute_and_save_rf_cells(datadir):
  runs = fam.get_successful_runs(datadir)
  print("Reading trees...")
  trees = read_all_trees(datadir, runs)
  print("Checking invalid runs...")
  #invalid_runs = get_invalid_runs(trees, runs)
  print("Removing invalid runs... (NOT IMPLEMENTED)")
  pass
  runs_to_compare = get_runs_to_compare(runs)
  print("Computing RF distances...") 
  rf_cells = compute_all_rf_cells(trees, runs_to_compare)
  print("Saving cells")
  save_rf_cells(datadir, rf_cells) 


def get_run_keys_from_rf_cells(rf_cells):
  run_keys = []
  for family in rf_cells:
    for run_key in rf_cells[family]:
      run_keys.append(run_key)
    break
  return run_keys

def export_metric(datadir, metric_dict, metric_name, benched_run):
  printer = AlignedPrinter()
  saved_metrics.save_dico(datadir, metric_dict, metric_name)
  for run_key in metric_dict:
    split = run_key.split(" - ")
    suffix = ""
    if (benched_run == split[0] or benched_run == split[1]):
      suffix += "\t <-- "
    printer.add("- " + run_key + ":",  str(metric_dict[run_key]) + suffix)
  printer.sort_right_float()
  printer.display()
  print("")


def compute_and_export_metrics(datadir, benched_run):
  rf_cells = load_rf_cells(datadir)
  run_keys = get_run_keys_from_rf_cells(rf_cells)
  total_nodes = {}
  total_rf = {}
  total_rrf = {}
  families_number = len(rf_cells)
  for run_key in run_keys:
    total_nodes[run_key] = 0.0
    total_rf[run_key] = 0.0
    total_rrf[run_key] = 0.0
  for family in rf_cells:
    family_rf_cells = rf_cells[family]
    for key in family_rf_cells:
      total_rf[key] += family_rf_cells[key][0]
      total_nodes[key] += family_rf_cells[key][1]
      total_rrf[key] += (family_rf_cells[key][0] / family_rf_cells[key][1])
 
  
  average_rrf = {}
  relative_arf = {}
  for key in run_keys:
    assert(total_nodes[key] > 0.0)
    average_rrf[key] = total_rrf[key] / families_number
    relative_arf[key] = total_rf[key] / total_nodes[key]
  print("Average relative RF:")
  export_metric(datadir, average_rrf, "average_rrf", benched_run)
  
  #print("Relative average RF:")
  #export_metric(datadir, relative_arf, "relative_arf")

def analyze(datadir, benched_run = "lastRun"):
  compute_and_save_rf_cells(datadir)
  compute_and_export_metrics(datadir, benched_run)



if __name__ == '__main__':
  if (len(sys.argv) < 2):
    print("Syntax: families_dir [benched_run]")
    exit(1)
  print(" ".join(sys.argv))
  datadir = sys.argv[1]
  benched_run = ""
  if (len(sys.argv) > 2):
    benched_run = sys.argv[2]
  analyze(datadir, benched_run)

