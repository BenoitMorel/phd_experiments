import sys
from ete3 import Tree


def contract_rec(node, min_bl):
  children = node.get_children()
  if (len(children) == 0):
    return
  if (not node.is_root() and node.dist <= min_bl):
    print("Contracting!")
    parent = node.get_ancestors()[0]
    assert(node in parent.get_children())
    for child in children:
      child.detach()
      parent.add_child(child)
    node.detach()
  for child in children:
    contract_rec(child, min_bl)
    

def contract_short_bl(input_tree, output_tree, min_bl):
  tree = Tree(input_tree)
  contract_rec(tree, min_bl)  
  tree.write(outfile = output_tree)


if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax: input_tree output_tree")
    exit(1)
  contract_short_bl(sys.argv[1], sys.argv[2], float(sys.argv[3]))





