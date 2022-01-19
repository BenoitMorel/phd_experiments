import sys
from ete3 import Tree


def contract_rec(node, min_bl, min_support):
  children = node.get_children()
  if (len(children) == 0):
    return
  bl = node.dist
  support = node.support
  is_polytomy = False
  is_polytomy  = node.dist <= min_bl or node.support < min_support
  is_polytomy = is_polytomy and not node.is_root()
  if (is_polytomy):
    parent = node.get_ancestors()[0]
    assert(node in parent.get_children())
    for child in children:
      child.detach()
      parent.add_child(child)
    node.detach()
  for child in children:
    contract_rec(child, min_bl, min_support)
    

def contract_branches(input_tree, output_tree, min_bl, min_support):
  tree = Tree(input_tree, format = 0)
  contract_rec(tree, min_bl, min_support)  
  tree.write(outfile = output_tree)


if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax: input_tree output_tree minbl minsupport")
    exit(1)
  contract_branches(sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]))





