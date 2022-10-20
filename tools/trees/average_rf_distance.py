import os
import sys
import read_tree
import rf_distance

def average_rf_distance(tree_path):
  trees = read_tree.read_trees_list(tree_path)
  s = 0.0
  for tree1 in trees:
    for tree2 in trees:
      d = rf_distance.get_relative_rf(tree1, tree2)
      s += d
  s /= float(len(trees) * len(trees))
  return s

if (__name__ == "__main__"):
  if (len(sys.argv) < 2):
    print("Syntax python " + os.path.basename(__file__) + " trees")
    sys.exit(1)
  trees = sys.argv[1]
  print(average_rf_distance(trees))
