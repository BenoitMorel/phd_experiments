from ete3 import Tree
import sys
import os


def cut_long_branches(tree_file, max_dist):
  tree = Tree(tree_file)
  nodes_to_cut = []
  trees = [tree]
  for node in tree.traverse():
    if (node.dist > max_dist):
      nodes_to_cut.append(node)
  for node in nodes_to_cut:
    trees.append(node.detach())
  return trees



if (__name__ == "__main__"):
  if (len(sys.argv) != 2): 
    print("Syntax: python " + os.path.basename(__file__) + " tree_file") 
    exit(1)

  tree_file = sys.argv[1]
  for tree in cut_long_branches(tree_file, 2.0):
    print(tree)

