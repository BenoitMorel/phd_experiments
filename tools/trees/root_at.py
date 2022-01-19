import re
import sys
import ete3

def root_at(input_tree, output_tree, taxon):
  tree = ete3.Tree(input_tree)
  leaf = None
  for node in tree.get_leaves():
    if (node.name == taxon):
      leaf = node
  assert(leaf != None)
  tree.set_outgroup(leaf)
  tree.write(outfile = output_tree)
  

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax: input_tree output_tree label")
    exit(1)
  root_at(sys.argv[1], sys.argv[2], sys.argv[3])




