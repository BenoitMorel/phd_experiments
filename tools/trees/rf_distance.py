import sys
import os
from ete3 import Tree

def str_2(ll):
  return "{0:.2f}".format(ll)

def str_4(ll):
  return "{0:.4f}".format(ll)

def ete3_rf(tree1, tree2):
  return tree1.robinson_foulds(tree2, unrooted_trees=True, skip_large_polytomies = True)

def get_rf(tree1, tree2):
  return ete3_rf(tree1, tree2)[0]

def get_relative_rf(tree1, tree2):
  rf = ete3_rf(tree1, tree2)
  return float(rf[0]) / float(rf[1])

def get_rf_distance_tuple(tree1_file, tree2_file):
  tree1 = Tree(tree1_file, format=1)
  tree2 = Tree(tree2_file, format=1)
  rf = ete3_rf(tree1, tree2)[0]
  return (float(rf[0]) / float(rf[1]), rf[0])

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python rf_distance.py tree1 tree2")
    sys.exit(1)
  tree1 = Tree(sys.argv[1], format=1)
  tree2 = Tree(sys.argv[2], format=1)
  print("Relative RF: " + str_4(get_relative_rf(tree1, tree2)))
  print("RF: " + str(get_rf(tree1, tree2)))
