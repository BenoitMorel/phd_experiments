
import sys
import os
from ete3 import Tree


def compute_induced_tree(tree_to_prune, subtree, output_tree):
  tree1 = Tree(tree_to_prune)
  leaf_set = Tree(subtree).get_leaf_names()
  tree1.prune(leaf_set)
  tree1.write(format = 1, outfile = output_tree)

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python diff.py tree_to_prune subtree output_tree")
    sys.exit(1)
  tree_to_prune = sys.argv[1]
  subtree = sys.argv[2]
  output_tree = sys.argv[3]
  compute_induced_tree(tree_to_prune, subtree, output_tree)





