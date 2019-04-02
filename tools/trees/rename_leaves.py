import sys
from ete3 import Tree

def rename_leaves(input_tree_str, dictionary):
  tree = Tree(input_tree_str, format=1)
  for leaf in tree.get_leaves():
    leaf.name = dictionary[leaf.name]
  return tree.write()

  

