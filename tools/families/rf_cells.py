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

def get_gene_trees_list(method, dataset_dir, msa):
  tree_path = fam.get_gene_tree(fam.get_family_path(dataset_dir, msa), method)
  return read_trees_list(tree_path)

def read_trees_for_family(datadir, family, methods):
  family_trees = {}
  for method in methods:
      family_trees[method] = get_gene_trees_list(method, datadir, family)
  return family_trees

def read_all_trees(datadir, methods):
  trees = {}
  families = fam.get_families_list(datadir)
  for family in families:
    trees[family] = read_trees_for_family(datadir, family, methods)
  return trees

def get_invalid_methods(trees, methods):
  invalid_methods = []
  for method in methods:
    for family in trees:
      if (not method in trees[family]):
        print("invalid method: " + method)
        invalid_methods.append(method)
        break
  return invalid_methods

def get_methods_to_compare(methods):
  methods_to_compare = []
  for method in methods:
    methods_to_compare.append((method, "true"))
  return methods_to_compare

def get_method_key(m1, m2):
  return m1 + " - " + m2

def compute_all_rf_cells(trees, methods_to_compare):
  rf_cells = {}
  for family in trees:
    family_trees = trees[family]
    rf_cells[family] = {}
    for method1, method2 in methods_to_compare:
      key = get_method_key(method1, method2)
      trees1 = family_trees[method1]
      trees2 = family_trees[method2]
      rf_cells[family][key] = rf_distance.ete3_average_rf_from_list(trees1, trees2) 
  return rf_cells

def save_rf_cells(datadir, rf_cells):
  output = fam.get_raw_rf_cells_file(datadir)
  pickle.dump(rf_cells, open(output, "wb"))

def load_rf_cells(datadir):
  return pickle.load(open(fam.get_raw_rf_cells_file(datadir), "rb"))

def compute_and_save_rf_cells(datadir):
  methods = fam.get_ran_methods(datadir)
  print("Reading trees...")
  trees = read_all_trees(datadir, methods)
  print("Checking invalid methods...")
  invalid_methods = get_invalid_methods(trees, methods)
  print("Removing invalid methods... (NOT IMPLEMENTED)")
  pass
  methods_to_compare = get_methods_to_compare(methods)
  print("Computing RF distances...") 
  rf_cells = compute_all_rf_cells(trees, methods_to_compare)
  print("Saving cells")
  save_rf_cells(datadir, rf_cells) 


def get_method_keys_from_rf_cells(rf_cells):
  method_keys = []
  for family in rf_cells:
    for method_key in rf_cells[family]:
      method_keys.append(method_key)
    break
  return method_keys

def export_metric(datadir, metric_dict, metric_name):
  printer = AlignedPrinter()
  saved_metrics.save_dico(datadir, metric_dict, metric_name)
  for method_key in metric_dict:
    printer.add("- " + method_key + ":",  str(metric_dict[method_key]))
  printer.sort_right_float()
  printer.display()
  print("")


def compute_and_export_metrics(datadir):
  rf_cells = load_rf_cells(datadir)
  method_keys = get_method_keys_from_rf_cells(rf_cells)
  total_nodes = {}
  total_rf = {}
  total_rrf = {}
  families_number = len(rf_cells)
  for method_key in method_keys:
    total_nodes[method_key] = 0.0
    total_rf[method_key] = 0.0
    total_rrf[method_key] = 0.0
  for family in rf_cells:
    family_rf_cells = rf_cells[family]
    for key in family_rf_cells:
      total_rf[key] += family_rf_cells[key][0]
      total_nodes[key] += family_rf_cells[key][1]
      total_rrf[key] += (family_rf_cells[key][0] / family_rf_cells[key][1])
 
  
  average_rrf = {}
  relative_arf = {}
  for key in method_keys:
    assert(total_nodes[key] > 0.0)
    average_rrf[key] = total_rrf[key] / families_number
    relative_arf[key] = total_rf[key] / total_nodes[key]
  print("Average relative RF:")
  export_metric(datadir, average_rrf, "average_rrf")
  print("Relative average RF:")
  export_metric(datadir, relative_arf, "relative_arf")




if __name__ == '__main__':
  if (len(sys.argv) < 2):
    print("Syntax: families_dir [method]")
    exit(1)
  print(" ".join(sys.argv))
  datadir = sys.argv[1]
  compute_and_save_rf_cells(datadir)
  compute_and_export_metrics(datadir)


