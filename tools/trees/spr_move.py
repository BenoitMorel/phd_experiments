
import sys
import os
import ete3
import read_tree

"""
  Applies an SPR move that prunes the node with label prune_label
  and regrafts it at the branch parent of the node with label
  regraft_label
  Only works if the pruned node and its parent are not the root
"""
def apply_spr(tree, prune_label, regraft_label):
  prune = None
  regraft = None
  for node in tree.traverse():
    if (node.name == prune_label):
      prune = node
    if (node.name == regraft_label):
      regraft = node
  if (prune == None):
    print("Can't find a prune node with label " + prune_label)
    return False
  if (regraft == None):
    print("Can't find a regraft node with label " + regraft_label)
    return False
  if (prune == tree or prune.up == tree):
    print("The prune node and its parent should not be the root of the tree")
    return False
  # we will also move the parent of prune
  parent_prune = prune.up
  grandparent_prune = parent_prune.up
  assert(len(prune.get_sisters()) == 1)
  sister_prune = prune.get_sisters()[0]
  parent_regraft = regraft.up
# do the pruning
  parent_prune = parent_prune.detach()
  sister_prune = sister_prune.detach()
  grandparent_prune.add_child(sister_prune)
# do the regraft
  regraft = regraft.detach()
  parent_regraft.add_child(parent_prune)
  parent_prune.add_child(regraft)



def generate_spr(input_tree, output_tree, prune_label, regraft_label):
  tree = read_tree.read_tree(input_tree)
  apply_spr(tree, prune_label, regraft_label)
  tree.write(outfile = output_tree)


if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax: input output prune_label regraft_label")
    exit(1)
  generate_spr(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


