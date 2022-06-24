import sys
import os
import ete3

def fast_prune_rec(node, to_keep):
  if (node.is_leaf()):
    if (node.name in to_keep):
      return node.name
    else:
      return ""
  children_strings = []
  for child in node.get_children():
    s = fast_prune_rec(child, to_keep)
    if (len(s) > 0):
      children_strings.append(s)
  if (len(children_strings) == 0):
    return ""
  if (len(children_strings) == 1):
    return children_strings[0]
  return "(" + ",".join(children_strings) + ")"


def fast_prune(tree, to_keep):
  s = fast_prune_rec(tree, to_keep) + ";"
  return ete3.Tree(s, format = 1)




