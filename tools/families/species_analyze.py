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
    if (len(split) == 4):
      res.append(split[0] + "." + split[1])
  return res

def get_species_tree(datadir, run):
  return os.path.join(fam.get_species_dir(datadir), run + ".speciesTree.newick")


def analyze(datadir):
  runs = get_runs(datadir)
  true_tree = read_tree(fam.get_species_tree(datadir))
  trees = {}
  for run in runs:
    trees[run] = read_tree(get_species_tree(datadir, run))

  print("")
  print("Rooted average RF:")
  rooted_printer = AlignedPrinter()
  for run in trees:
    try:
      tree = trees[run]
      rooted_rf_cell = true_tree.robinson_foulds(tree, unrooted_trees=False, correct_by_polytomy_size = True)
      rooted_arf = float(rooted_rf_cell[0]) / float(rooted_rf_cell[1])
      if (len(tree) != len(true_tree)):
        rooted_arf = 10000000.0
      saved_metrics.save_metrics(datadir, run, str(rooted_arf), "species_rooted_rf") 
      rooted_printer.add(run + ":", str(rooted_arf))
    except:
      pass
  rooted_printer.sort_right_float()
  rooted_printer.display()

  print("")
  print("Unrooted average RF:")
  unrooted_printer = AlignedPrinter()
  for run in trees:
    tree = trees[run]
    if (tree == None):
      continue
    unrooted_rf_cell = true_tree.robinson_foulds(tree, unrooted_trees=True, correct_by_polytomy_size = True)
    unrooted_arf = float(unrooted_rf_cell[0]) / float(unrooted_rf_cell[1])
    if (len(tree) != len(true_tree)):
      unrooted_arf = 10000000.0
    saved_metrics.save_metrics(datadir, run, str(unrooted_arf), "species_unrooted_rf") 
    unrooted_printer.add(run + ":", str(unrooted_arf))
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


