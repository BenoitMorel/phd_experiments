import sys
import os
from ete3 import Tree
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import experiments as exp
import fam


def get_dict(datadir):
  d = {}
  for line in open(fam.get_species_dict(datadir)).readlines():
    sp = line.replace("\n", "").split(":")
    d[sp[0]] = sp[1]
  return  d

def translate(species_tree_path):
  tree = Tree(species_tree_path, format = 1)
  datadir = species_tree_path.split("species_tree")[0]
  species_dict = get_dict(datadir)
  for leaf in tree.get_leaves():
    leaf.name = species_dict[leaf.name]
  print(tree.write())


if (__name__ == "__main__"): 
  if (len(sys.argv) < 2): 
    print("Syntax: python " + os.path.basename(__file__) + " species_tree")
    exit(1)

translate(sys.argv[1])


