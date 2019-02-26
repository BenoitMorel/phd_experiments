import os
import sys
from ete3 import Tree
from read_tree import read_tree
import random

def _resolve(node):
  if len(node.children) > 2:
    children = list(node.children)
    random.shuffle(children)    
    node.children = []
    next_node = root = node
    for i in range(len(children)-2):
      next_node = next_node.add_child()
      next_node.dist = 1.0
      next_node.support = 0.0

    next_node = root
    for ch in children:
      next_node.add_child(ch)
      if ch != children[-2]:
         next_node = next_node.children[0]

def resolve_polytomies(tree):
  target = [tree]
  target.extend([n for n in tree.get_descendants()])
  for n in target:
    _resolve(n)


def make_binary(input_tree, output_tree, seed):
  random.seed(seed)
  tree = read_tree(input_tree)
  tree.resolve_polytomy()
  tree = read_tree(input_tree)
  resolve_polytomies(tree)
  with open(output_tree, "w") as writer:
    tree.write(outfile = output_tree)

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax: input_tree output_tree seed")
    exit(1)
  make_binary(sys.argv[1], sys.argv[2], int(sys.argv[3]))

