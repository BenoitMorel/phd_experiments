import re
import sys
from ete3 import Tree

def read_tree(tree):
  lines = open(tree).readlines()
  for line in lines:
    if (not line.startswith(">")):
      return Tree(line, format=1)
  return None

def read_trees(trees):
  lines = open(trees).readlines()
  res = []
  for line in lines:
    if (not line.startswith(">")):
      res.append(Tree(line, format=1))
  return res

def save_trees(trees, output_files):
  with open(output_files, "w") as writer:
    for tree in trees:
      writer.write(tree.write() + "\n")


def cut_keep_first_elems(input_tree, output_tree, separator, keep_number):
  tree = read_tree(input_tree)
  for leaf in tree.get_leaves():
    leaf.name = separator.join(leaf.name.split(separator)[0:keep_number])
  tree.write(outfile = output_tree) 

def remove_prefix_from_trees(input_trees, output_trees, separator = "_"):
  trees = read_trees(input_trees)
  for tree in trees:
    for leaf in tree.get_leaves():
      leaf.name = separator.join(leaf.name.split(separator)[1:])
  save_trees(trees, output_trees)


#if (__name__ == "__main__"):
#  if (len(sys.argv) != 3):
#    print("Syntax: input_tree output_tree")
#    exit(1)
#  cut_node_names(sys.argv[1], sys.argv[2])




