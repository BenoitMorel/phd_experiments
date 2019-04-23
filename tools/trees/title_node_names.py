import re
import sys
from ete3 import Tree

def read_tree(tree):
  lines = open(tree).readlines()
  for line in lines:
    if (not line.startswith(">")):
      return Tree(line, format=1)
  return None

def get_title(name):
  split = name.split("_")
  for i in range(0, len(split)):
    split[i] = split[i].title()
  return "UUU".join(split)

def title_node_names(input_tree, output_tree):
  tree = read_tree(input_tree)
  for leaf in tree.get_leaves():
    leaf.name = get_title(leaf.name)
  tree.write(outfile = output_tree) 

if (__name__ == "__main__"): 
  if (len(sys.argv) != 3): 
     print("Syntax: python " + os.path.basename(__file__) + " input_tree output_tree")
     exit(1)
  title_node_names(sys.argv[1], sys.argv[2])


