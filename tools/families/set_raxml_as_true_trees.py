import os
import sys
import shutil

if (len(sys.argv) != 2):
  print("Syntax: python set_raxml_as_true_tree.py dataset_dir")

dataset_dir = sys.argv[1]

families_dir = os.path.join(dataset_dir, "families")
for family in os.listdir(families_dir):
  raxml_tree = os.path.join(families_dir, family, "raxmlGeneTree.newick")
  true_tree = os.path.join(families_dir, family, "trueGeneTree.newick")
  shutil.copy(raxml_tree, true_tree)

