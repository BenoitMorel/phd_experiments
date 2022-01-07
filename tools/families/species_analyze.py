import sys
import os
from ete3 import Tree
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "print"))
sys.path.insert(0, os.path.join("script"))
from read_tree import read_tree
import saved_metrics
import fam
from aligned_printer import AlignedPrinter
import experiments as exp

def get_runs(datadir):
  res = []
  for f in os.listdir(fam.get_species_dir(datadir)):
    split = f.split(".")
    run = ".".join(split[0:-3])
    if (len(split) >= 4):
      res.append(run + "." + split[-3])
  return res

def get_species_tree(datadir, run):
  return os.path.join(fam.get_species_dir(datadir), run + ".speciesTree.newick")

def get_split(tree):
  res = []
  for child in tree.get_children():
    res.append(set(child.get_leaf_names()))
  return res

def get_split_score(sp1, sp2):
  if (len(sp1) != 2):
    return 1.0
  if (len(sp2) != 2):
    return 1.0
  if (sp1[0] == sp2[0]):
    return 0.0
  if (sp1[0] == sp2[1]):
    return 0.0
  return 1.0

def analyze(datadir):
  runs = get_runs(datadir)
  true_tree = read_tree(fam.get_species_tree(datadir))
  true_split = get_split(true_tree)
  trees = {}
  for run in runs:
    try:
      trees[run] = read_tree(get_species_tree(datadir, run))
    except:
      print("Cannot read " + get_species_tree(datadir, run))
      #raise

  root_distance_printer = AlignedPrinter()
  rooted_printer = AlignedPrinter()
  unrooted_printer = AlignedPrinter()

  for run in trees:
    #try:
      tree = trees[run]
      if (tree == None):
        continue
      # add a fake root if the tree is unrooted
      if (len(tree.children) == 3):
        tree.set_outgroup(tree.children[0])
      # check that split is correct
      tree_split = get_split(tree)
      split_score = get_split_score(true_split, tree_split)
      saved_metrics.save_metrics(datadir, run, str(split_score), "root_split") 

      # compute RF distances
      rooted_rf = [0.0,0.0]
      try:
        rooted_rf = true_tree.robinson_foulds(tree, unrooted_trees=False, correct_by_polytomy_size = True)
      except:
        pass
      unrooted_rf = true_tree.robinson_foulds(tree, unrooted_trees= True, correct_by_polytomy_size = True)
      # relativee distances 
      rooted_arf = 10000000.0 
      unrooted_arf = 10000000.0 
      root_distance = 10000000.0
      if (len(tree) == len(true_tree)):
        if (float(rooted_rf[1] != 0.0)):
          rooted_arf = float(rooted_rf[0]) / float(rooted_rf[1])
          root_distance = float(rooted_rf[0] - unrooted_rf[0]) / float(rooted_rf[1])
          
        if (float(unrooted_rf[1] != 0.0)):
          unrooted_arf = float(unrooted_rf[0]) / float(unrooted_rf[1])
      else:
        print("ERROR " + str(len(tree)) + " " + str(len(true_tree)))
      # save metrics
      saved_metrics.save_metrics(datadir, run, str(rooted_arf), "species_rooted_rf") 
      saved_metrics.save_metrics(datadir, run, str(unrooted_arf), "species_unrooted_rf") 
      saved_metrics.save_metrics(datadir, run, str(root_distance), "root_distance") 
      # add to printer
      rooted_printer.add(run + ":", str(rooted_arf))
      unrooted_printer.add(run + ":", str(unrooted_arf))
      root_distance_printer.add(run + ":" , str(root_distance))
    #except:
   #   pass

  # print results
  print("")
  print("Root distance:")
  root_distance_printer.sort_right_float()
  root_distance_printer.display()

  print("")
  print("Rooted average RF:")
  rooted_printer.sort_right_float()
  rooted_printer.display()

  print("")
  print("Unrooted average RF:")
  unrooted_printer.sort_right_float()
  unrooted_printer.display()
    

def analyze_all():
  for datadir in os.listdir(exp.families_datasets_root):
    datadir = os.path.join(exp.families_datasets_root, datadir)
    try:
      analyze(datadir)
      print("Analyzed " + datadir)
    except:
      pass

if __name__ == '__main__':
  if (len(sys.argv) < 2):
    print("Syntax: families_dir")
    exit(1)
  print(" ".join(sys.argv))
  datadir = sys.argv[1]
  if (datadir == "all"):
    analyze_all()
  else:
    analyze(datadir)


