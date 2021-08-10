import sys
import os
from ete3 import Tree
import ete3
print("hey  " + ete3.__file__)

def read_tree(tree):
  lines = open(tree).readlines()
  for line in lines:
    if (not line.startswith(">")):
      return Tree(line, format=1)
  return None

def rescale_bl(input_tree, output_tree, scale):
  tree = read_tree(input_tree)
  for node in tree.traverse("postorder"):
    node.dist *= scale
  tree.write(outfile=output_tree, format = 5, dist_formatter= "%.15f")

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax: python rescale_bl.py input_tree output_tree scale")
    exit(1)
  input_tree = sys.argv[1]
  output_tree = sys.argv[2]
  scale = float(sys.argv[3])
  rescale_bl(input_tree, output_tree, scale)
