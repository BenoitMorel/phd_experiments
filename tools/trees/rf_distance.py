import sys
import os
from ete3 import Tree



def str_2(ll):
  return "{0:.2f}".format(ll)

def str_4(ll):
  return "{0:.4f}".format(ll)

def ete3_rf(tree1, tree2):
  return tree1.robinson_foulds(tree2, unrooted_trees=False, skip_large_polytomies = True)

def ete3_average_rf_from_list(tree_list_1, tree_list_2, rooted):
  average_cell = [0.0, 0.0]
  for tree1 in tree_list_1:
    for tree2 in tree_list_2:
      cell = tree1.robinson_foulds(tree2, unrooted_trees= (not rooted))
      average_cell[0] += float(cell[0])
      average_cell[1] += float(cell[1])
  average_cell[0] /= float(len(tree_list_1) * len(tree_list_2))
  average_cell[1] /= float(len(tree_list_1) * len(tree_list_2))
  return average_cell

def ete3_min_rf_from_list(tree_list_1, tree_list_2):
  min_cell = [1000000000, 110000000]
  for tree1 in tree_list_1:
    for tree2 in tree_list_2:
      cell = tree1.robinson_foulds(tree2, unrooted_trees=True)
      if (cell[0] < min_cell[0]):
        min_cell = cell
  return min_cell

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
