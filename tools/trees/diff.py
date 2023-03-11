import sys
import os
import ete3

def tree_diff(tree_file_1, tree_file_2):
  tree1 = ete3.Tree(tree_file_1, format = 1)
  tree2 = ete3.Tree(tree_file_2, format = 1)
  leaves1 = set(tree1.get_leaf_names())
  leaves2 = set(tree2.get_leaf_names())
  only1 = leaves1 - leaves2
  only2 = leaves2 - leaves1
  print("")
  print("Leaves that are in first tree only:")
  print(" ".join(sorted(list(only1))))
  print("")
  print("")
  print("Leaves that are in second tree only:")
  print(" ".join(sorted(list(only2))))
  print("")

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python diff.py tree1 tree2")
    sys.exit(1)
  tree_file_1 = sys.argv[1]
  tree_file_2 = sys.argv[2]
  tree_diff(tree_file_1, tree_file_2)




