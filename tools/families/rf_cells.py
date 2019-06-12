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
  print("Cells: " + str(rf_cells))
  save_rf_cells(datadir, rf_cells) 

if __name__ == '__main__':
  if (len(sys.argv) < 2):
    print("Syntax: families_dir [method]")
    exit(1)
  print(" ".join(sys.argv))
  dataset_dir = sys.argv[1]
  compute_and_save_rf_cells(dataset_dir)



