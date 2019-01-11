import sys
import os
from ete3 import Tree

def read_tree(tree):
  lines = open(tree).readlines()
  for line in lines:
    if (not line.startswith(">")):
      return Tree(line, format=1)
  return None

def get_rf(tree1, tree2):
  return tree1.robinson_foulds(tree2, unrooted_trees=True)[0]

def get_relative_rf(tree1, tree2):
  rf = tree1.robinson_foulds(tree2, unrooted_trees=True)
  return float(rf[0]) / float(rf[1])

def get_nodes(tree1):
  rf = tree1.robinson_foulds(tree1, unrooted_trees=True)
  return rf[1]



def analyse(dataset_dir, pargenes_dir):
  analysed_msas = 0
  total_nodes_number = 0
  methods = ["Raxml-ng", "Treerecs", "Phyldog", "JointSearch"]
  methods_tree_files = {}
  methods_tree_files["Raxml-ng"] = "raxmlGeneTree.newick"
  methods_tree_files["Treerecs"] = "treerecsGeneTree.newick"
  methods_tree_files["Phyldog"] = "phyldogGeneTree.newick"
  methods_tree_files["JointSearch"] = "jointsearch.newick"
  methods_trees_number = {}
  for m in methods:
    methods_trees_number[m] = 0
  total_rrf = {}
  total_rf = {}
  true_matches = {}
  best_tree = {}
  for method in methods:
    total_rrf[method] = 0.0
    total_rf[method] = 0.0
    true_matches[method] = 0
    best_tree[method] = 0
  for msa in os.listdir(dataset_dir):   
    trees = {}
    family_path = os.path.join(dataset_dir, msa)
    try:
      true_tree = Tree(os.path.join(family_path, "trueGeneTree.newick"), format=1) 
    except:
      continue
    ok = True
    for method in methods:
      if (method == "JointSearch"):
        prefix = os.path.join(pargenes_dir, "results", msa)
      else:
        prefix = family_path
      try:
        trees[method] = read_tree(os.path.join(prefix, methods_tree_files[method]))
        methods_trees_number[method] += 1
      except:
        try:
          #trees[method] = read_tree(os.path.join(prefix, methods_tree_files["Raxml-ng"]))
          trees[method] = Tree(os.path.join(family_path, "trueGeneTree.newick"), format=1)
        except:
          trees[method] = None
    if (not ok):
      continue
    #newmethods = []
    #for m in methods:
    #  if (methods_trees_number[m] == 0):
    #    print("Remove method " + m  + " from analysis")
    #  else:
    #    newmethods.append(m)
    #methods = newmethods
    best_rrf = 1
    rrf = {}
    rf = {}
    analysed_msas += 1
    for method in methods:
      if (trees[method] == None):
        continue
      rf_cell = trees[method].robinson_foulds(true_tree, unrooted_trees=True)
      rrf[method] = float(rf_cell[0]) / float(rf_cell[1])
      rf[method] = float(rf_cell[0])
      best_rrf = min(best_rrf, rrf[method])
      total_rrf[method] += rrf[method]
      total_rf[method] += float(rf_cell[0])
      if (rf_cell[0] == 0):
        true_matches[method] += 1
      if (method == methods[0]): #do it only once
        total_nodes_number += rf_cell[1]
    for method in methods:
      if (best_rrf == rrf[method]):
        best_tree[method] += 1

  if (analysed_msas == 0):
    print("did not manage to analyse any MSA")
    exit(1)
  print("Number of gene families: " + str(analysed_msas))
  print("")
  print("Average (over the gene families) relative RF distance to the true trees:")
  for method in methods:
    print("- " + method + ":\t" + str(total_rrf[method] / float(analysed_msas)))
  print("")
  
  print("Normalized average (over the gene families) RF distance to the true trees:")
  for method in methods:
    print("- " + method + ":\t" + str(total_rf[method] / float(total_nodes_number)))
  print("")
  
  print("Number of gene families for which a method reaches the smallest relative RF to the true trees compared with the other methods:")
  for method in methods:
    print("- " + method + ":\t" + str(best_tree[method]) + "/" + str(analysed_msas))
  print("")
  
  print("Number of gene families for which a method finds the true tree:")
  for method in methods:
    print("- " + method + ":\t" + str(true_matches[method]) + "/" + str(analysed_msas))
  print("")

  print("Analysed trees:")
  for method in methods:
    print("- " + method + ":\t" + str(methods_trees_number[method]))

if __name__ == '__main__':
  if (len(sys.argv) != 3):
    print("Syntax: families_dir pargenes_dir")
    exit(1)
  dataset_dir = sys.argv[1]
  pargenes_dir = sys.argv[2]
  analyse(dataset_dir, pargenes_dir)



