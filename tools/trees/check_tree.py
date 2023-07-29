import os
import sys
from ete3 import Tree
from read_tree import read_tree


def check(tree_file):
  tree = read_tree(tree_file)
  leaves = set()
  for node in tree.get_leaves():
    if (node.name in leaves):
      print("DUPLICATE " + node.name)
    leaves.add(node.name) 


if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax: input_tree")
    exit(1)
  check(sys.argv[1])

