import sys
import os

from ete3 import Tree

def str_2(ll):
  return "{0:.2f}".format(ll)

def str_4(ll):
  return "{0:.4f}".format(ll)

def get_rf(tree1, tree2):
  return tree1.robinson_foulds(tree2, unrooted_trees=True)[0]

def get_relative_rf(tree1, tree2):
  rf = tree1.robinson_foulds(tree2, unrooted_trees=True)
  return float(rf[0]) / float(rf[1])

if (len(sys.argv) != 3):
  print("Syntax python rf_distance.py tree1 tree2")
  sys.exit(1)

tree1 = Tree(sys.argv[1], format=1)
tree2 = Tree(sys.argv[2], format=1)

print("Relative RF: " + str_4(get_relative_rf(tree1, tree2)))
print("RF: " + str(get_rf(tree1, tree2)))
