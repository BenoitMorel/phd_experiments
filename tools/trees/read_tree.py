import sys
import os
from ete3 import Tree

def read_tree(tree):
  with open(tree) as reader:
    while (True):
      line = reader.readline()
      if (not line):
        return None
      if (not line.startswith(">")):
        return Tree(line, format=1)
  return None


def read_trees_list(newick_file):
  trees = []
  lines = open(newick_file).readlines()
  for line in lines:
    if (not line.startswith(">")):
      trees.append(Tree(line, format=1))
  return trees


