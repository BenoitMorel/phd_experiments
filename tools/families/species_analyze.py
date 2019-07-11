import sys
import os
from ete3 import Tree
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "families"))
from read_tree import read_tree
import saved_metrics
import fam


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

  print("Rooted average RF:")
  for run in trees:
    try:
      tree = trees[run]
      rooted_rf_cell = true_tree.robinson_foulds(tree, unrooted_trees=False)
      rooted_arf = float(rooted_rf_cell[0]) / float(rooted_rf_cell[1])
      print(run + ":\t" + str(rooted_arf)) 
    except:
      pass

  print("Unrooted average RF:")
  for run in trees:
    tree = trees[run]
    unrooted_rf_cell = true_tree.robinson_foulds(tree, unrooted_trees=True)
    unrooted_arf = float(unrooted_rf_cell[0]) / float(unrooted_rf_cell[1])
    print(run + ":\t" + str(unrooted_arf)) 
    



if __name__ == '__main__':
  if (len(sys.argv) < 2):
    print("Syntax: families_dir")
    exit(1)
  print(" ".join(sys.argv))
  datadir = sys.argv[1]
  analyze(datadir)


