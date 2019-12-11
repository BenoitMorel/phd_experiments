import sys
import ete3



def remove_leaves(input_tree, output_tree, leaves):
  tree = ete3.Tree(input_tree, format=1)
  all_leaves = tree.get_leaf_names()
  leaves_to_keep = list(set(all_leaves) - set(leaves))
  tree.prune(leaves_to_keep)
  open(output_tree, "w").write(tree.write())


if (__name__ == "__main__"):
  if (len(sys.argv) < 4):
    print("Syntax python remove_leaves input output leave1 leave2...")
    sys.exit(1)
  input_tree = sys.argv[1]
  output_tree = sys.argv[2]
  leaves = sys.argv[3:]
  remove_leaves(input_tree, output_tree, leaves)

