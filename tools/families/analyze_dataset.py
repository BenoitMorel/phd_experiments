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

def get_nodes(tree1):
  rf = ete3_rf(tree1, tree1)
  return rf[1]


methods_tree_files = {}
methods_tree_files["True"] = "trueGeneTree.newick"
methods_tree_files["RAxML-NG"] = "raxmlGeneTree.newick"
methods_tree_files["Treerecs"] = "treerecsGeneTree.newick"
methods_tree_files["Phyldog"] = "phyldogGeneTree.newick"
methods_tree_files["Notung80"] = "notung80GeneTree.newick"
methods_tree_files["Notung101"] = "notung101GeneTree.newick"
methods_tree_files["ALE-DL"] = "ALE-DLGeneTree.newick"
methods_tree_files["ALE-DTL"] = "ALE-DTLGeneTree.newick"

class AlignedPrinter:
  def __init__(self):
    self._lefts = []
    self._rights = []

  def add(self, left, right):
    self._lefts.append(left)
    self._rights.append(right)

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
  tree_path = ""
  if (method in methods_tree_files):
    if (method == "True"):
      tree_path = os.path.join(dataset_dir, msa, methods_tree_files[method])
    else:
      tree_path = os.path.join(dataset_dir, msa, "gene_trees", methods_tree_files[method])
  else:
    tree_path = os.path.join(dataset_dir, msa, "results", method + ".newick") 
  return read_trees_list(tree_path)

def add_ran_methods(methods, families_dir):
  runs_dir = os.path.join(families_dir, os.listdir(families_dir)[0], "results")
  if (not os.path.isdir(runs_dir)):
    return
  for method in os.listdir(runs_dir):
    method = method.split(".")[0]
    #if (method != "lastRun"):
    methods.append(method)

def analyze_msa(params):
  msa, dataset_dir, methods, methods_to_compare = params
  trees = {}
  family_path = os.path.join(dataset_dir, msa)
  invalid_methods = []
  best_tree = {}
  
  for method in methods:
    try: 
      trees[method] = get_gene_trees_list(method, dataset_dir, msa)
    except:
      invalid_methods.append(method)
  for method in invalid_methods:
    if (method in methods):
      methods.remove(method)
    methods_to_compare = [p for p in methods_to_compare if (p[0] != method and p[1] != method)]
  best_rrf = 1
  rrf = {}
  for method_pair in methods_to_compare:
    method1 = method_pair[0]
    method2 = method_pair[1]
    methods_key = method1 + " - " + method2
    rf_cell = rf_distance.ete3_average_rf_from_list(trees[method1], trees[method2])
    if (rf_cell[1] == 0):
      print("null cell for " + methods_key + " " + msa)
      print(trees[method1])
      print(trees[method2])
      exit(1)
    rrf[methods_key] = float(rf_cell[0]) / float(rf_cell[1])
    if (method1 == "True" or method2 == "True"):
      best_rrf = min(best_rrf, rrf[methods_key])
  for method_pair in methods_to_compare:
    method1 = method_pair[0]
    method2 = method_pair[1]
    methods_key = method1 + " - " + method2
    if (method1 == "True"):
      if (best_rrf == rrf[methods_key]):
        best_tree[method2] = 1
    elif (method2 == "True"):
      if (best_rrf == rrf[methods_key]):
        best_tree[method1] = 1
  return rrf, best_tree, invalid_methods

def analyze(dataset_dir, benched_method = ""):
  print("To re-run: python " + os.path.realpath(__file__) + " " + dataset_dir + " " + benched_method)
  families_dir = os.path.join(dataset_dir, "families") 
  analyzed_msas = 0
  total_nodes_number = 0
  methods = ["True", "RAxML-NG", "Phyldog","Treerecs", "Phyldog", "Notung80", "Notung101", "ALE-DL", "ALE-DTL"]
  add_ran_methods(methods, families_dir)
  methods_trees_number = {}
  methods_to_compare = []
  for method in methods:
    if (method == "True"):
      continue
    methods_to_compare.append(("True", method))
  for m in methods:
    methods_trees_number[m] = 0
  total_rrf = {}
  best_tree = {}
  for method in methods:
    best_tree[method] = 0
  for method_pair in methods_to_compare:
    method1 = method_pair[0]
    method2 = method_pair[1]
    methods_key = method1 + " - " + method2
    total_rrf[methods_key] = 0.0
  analyze_msa_params = []
  #with concurrent.futures.ProcessPoolExecutor(1) as executor:
  for msa in os.listdir(families_dir):  
    rrf, bt, invalid_methods = analyze_msa((msa, families_dir, methods, methods_to_compare))
    for method in invalid_methods:
      print("remove " + method)
      print(methods_to_compare)
      methods_to_compare = [p for p in methods_to_compare if (p[0] != method and p[1] != method)]
      print(methods_to_compare)
      #del best_tree[method]
    for m in rrf:
      total_rrf[m] += rrf[m]
    for m in bt:
      best_tree[m] += bt[m]
    analyzed_msas += 1
  if (analyzed_msas == 0):
    print("did not manage to analyze any MSA")
    exit(1)
        
  

  print("Number of gene families: " + str(analyzed_msas))
  print("")

  print("Average (over the gene families) relative RF distances:")
  rrf_printer = AlignedPrinter()
  for method_pair in methods_to_compare:
    method1 = method_pair[0]
    method2 = method_pair[1]
    methods_key = method1 + " - " + method2
    arrow = ""
    if (method1 == "True" and method2 == benched_method):
      arrow = " <-- "
    average_rrf = str(total_rrf[methods_key] / float(analyzed_msas))
    if (method1 == "True"):
      saved_metrics.save_metrics(dataset_dir, method2, average_rrf, "average_rrf")
    rrf_printer.add("- " + methods_key + ":",  average_rrf + arrow)
  rrf_printer.display()
  print("")
  
  
  print("Number of gene families for which a method reaches the smallest relative RF to the true trees compared with the other methods:")
  best_printer = AlignedPrinter()
  for method in methods:
    if (method == "True"):
      continue
    best_printer.add("- " + method + ":",  str(best_tree[method]) + "/" + str(analyzed_msas))
  best_printer.display()
  print("")
  


if __name__ == '__main__':
  if (len(sys.argv) < 2):
    print("Syntax: families_dir [method]")
    exit(1)
  print(" ".join(sys.argv))
  dataset_dir = sys.argv[1]
  method = ""
  if (len(sys.argv) > 2):
    method = sys.argv[2]
  analyze(dataset_dir, method)



