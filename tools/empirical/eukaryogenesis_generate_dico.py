import sys
import os
sys.path.insert(0, 'tools/trees')
from read_tree import read_tree


def gen_dico(tree_long, tree_short, output):
  leaves_long = set(read_tree(tree_long).get_leaf_names())
  leaves_short = set(read_tree(tree_short).get_leaf_names())
  with open(output, "w") as writer:
    for label in leaves_long:
      short = label.split("|")[-1]
      assert(short in  leaves_short)
      writer.write(short + ":" + label + "\n")

if (__name__ == "__main__"):
  if (len(sys.argv) < 4):
    print("Syntax python " + os.path.basename(__file__) + " tree_long tree_short output_dico")
    sys.exit(1)
  tree_long = sys.argv[1]
  tree_short = sys.argv[2]
  output = sys.argv[3]
  gen_dico(tree_long, tree_short, output)
