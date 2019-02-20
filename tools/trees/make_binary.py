import os
import sys
from ete3 import Tree

def read_tree(tree):
  lines = open(tree).readlines()
  for line in lines:
    if (not line.startswith(">")):
      return Tree(line, format=1)
  return None

def make_binary(input_tree, output_tree):
  tree = read_tree(input_tree)
  tree.resolve_polytomy()
  tree = read_tree(input_tree)
  tree.resolve_polytomy()
  with open(output_tree, "w") as writer:
    tree.write(outfile = output_tree)

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax: input_tree output_tree")
    exit(1)
  make_binary(sys.argv[1], sys.argv[2])

