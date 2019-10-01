import sys
import os
from ete3 import Tree
from read_tree import read_tree

def analyze_dimensions(tree):
  print("Taxa number: " + str(len(tree.get_leaves())))

def get_nodes_number(tree_file):
  return len(read_tree(tree_file).get_leaves())

def analyze_polytomies(tree):
  print("Polytomies:")
  for node in tree.traverse("postorder"):
    children = node.get_children()
    if (len(children) > 2):
      only_leaves = True
      for child in children:
        if (not child.is_leaf()):
          only_leaves = False
          break
      print("  Polytomy of size " + str(len(children)))
      if (only_leaves):
        print("    all the children in this polytomy are leaves")
      else:
        print("    some children in this polytomy are internal nodes")
      print("    number of taxa under the polytomy: " + str(len(node.get_leaves())))

def analyze_tree(tree_file):
  tree = read_tree(tree_file)
  analyze_dimensions(tree)
  analyze_polytomies(tree)
 
def get_tree_taxa_number(tree_file):
  return len(read_tree(tree_file).get_leaves())

def is_ultrametric(tree_file, epsilon = 0.001):
  tree = read_tree(tree_file)
  ref = -1.0
  for node in tree.traverse("postorder"):
    if (node.is_leaf()):
      if (ref == -1.0):
        ref = node.get_distance(tree)
      else:
        if (abs(ref - node.get_distance(tree)) > epsilon):
          return False
  return True

def check_ultrametric_and_get_length(tree_file, epsilon = 0.001):
  tree = read_tree(tree_file)
  ref = -1.0
  for node in tree.traverse("postorder"):
    if (node.is_leaf()):
      if (ref == -1.0):
        ref = node.get_distance(tree)
      else:
        assert(abs(ref - node.get_distance(tree)) < epsilon)
  return ref

if (__name__ == "__main__"):
  if (len(sys.argv) != 2): 
    print("Syntax: python " + os.path.basename(__file__) + " tree_file") 
    exit(1)

  tree_file = sys.argv[1]
  analyze_tree(tree_file)
  if (is_ultrametric(tree_file)):
    print("Ultrametric")
  else:
    print("Not ultrametric")

