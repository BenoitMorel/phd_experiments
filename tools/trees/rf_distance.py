import sys
import os
from ete3 import Tree
import dendropy


def str_2(ll):
  return "{0:.2f}".format(ll)

def str_4(ll):
  return "{0:.4f}".format(ll)

def ete3_rf(tree1, tree2, unrooted = True):
  if (len(tree2.children) == 3):
    tree2.set_outgroup(tree2.children[0])
  if (len(tree1.children) == 3):
    tree1.set_outgroup(tree1.children[0])

  return tree1.robinson_foulds(tree2, unrooted_trees=unrooted, skip_large_polytomies = True, correct_by_polytomy_size = True)

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

def get_relative_rf(tree1, tree2, unrooted = True):
  
  rf = ete3_rf(tree1, tree2, unrooted)
  #print(rf[0])
  return float(rf[0]) / float(rf[1])

def get_rf_distance_tuple(tree1_file, tree2_file):
  tree1 = Tree(tree1_file, format=1)
  tree2 = Tree(tree2_file, format=1)
  rf = ete3_rf(tree1, tree2)
  return (float(rf[0]) / float(rf[1]), rf[0])

def dendropy_rf_distance(tree1_file, tree2_file):
  dataSet = dendropy.DataSet()
  dataSet.read(file=open(tree1_file), schema = "newick")
  dataSet.read(file=open(tree2_file), schema = "newick", taxon_namespace=dataSet.tree_lists[0].taxon_namespace)
  distance = dendropy.calculate.treecompare.symmetric_difference(dataSet.tree_lists[0][0], dataSet.tree_lists[1][0], is_bipartitions_updated=False)

  return float(distance)

if (__name__ == "__main__"):
  if (len(sys.argv) < 3):
    print("Syntax python rf_distance.py tree1 tree2 [rooted]")
    sys.exit(1)
  tree1 = Tree(sys.argv[1], format=1)
  tree2 = Tree(sys.argv[2], format=1)
  rooted = False
  if (len(sys.argv) > 3):
    rooted = int(sys.argv[3]) > 0
  if (rooted):
    print("ROOTED")
  print("Relative RF: " + str_4(get_relative_rf(tree1, tree2, not rooted)))
  #print("RF: " + str(get_rf(tree1, tree2)))
#  print("dendrop distance: " + str(dendropy_rf_distance(sys.argv[1], sys.argv[2])))

