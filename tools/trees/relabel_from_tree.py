import os
import sys
import read_tree
sys.path.insert(0, 'tools/trees')


"""
  This function assumes that both input trees have the exact
  same rooted topology, and only differ by their branch lengths
  and internal node labels.

  It loads tree_to_relabel_path, changes its internal node
  labels to match the ones from ref_tree_path and returns the
  new tree
"""
def load_relabeled(tree_to_relabel_path, ref_tree_path):
  ref_tree = read_tree.read_tree(ref_tree_path)
  ref_label_to_node = {}
  for node in ref_tree.traverse():
    if (not node.is_root()):
      ref_label_to_node[node.name] = node

  tree = read_tree.read_tree(tree_to_relabel_path)
  for node in tree.traverse("postorder"):
    if (not node.is_root()):
      ref_node = ref_label_to_node[node.name]
      node.up.name = ref_node.up.name
  return tree



if (__name__ == "__main__"):
  if (len(sys.argv) < 4):
    print("Syntax python " + os.path.basename(__file__) + "tree_to_relabel_path, ref_tree_path, output")
    sys.exit(1)
  tree_to_relabel_path = sys.argv[1]
  ref_tree_path = sys.argv[2]
  output = sys.argv[3]
  tree = load_relabeled(tree_to_relabel_path, ref_tree_path)
  tree.write(outfile=output, format = 1) 
  

