import os
import sys
import shutil


def extract(dataset, method_tag):
  families_path = os.path.join(dataset, "families")
  stag_path = os.path.join(dataset, "stag", method_tag + "_gene_trees")
  try:
    os.makedirs(stag_path)
  except:
    pass
  for family in os.listdir(families_path):
    gene_tree = os.path.join(families_path, family, method_tag + "GeneTree.newick")
    dest = os.path.join(stag_path, family + ".newick")
    shutil.copyfile(gene_tree, dest)
  return stag_path

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax: python extract_gene_trees.py dataset method_tag")
    exit(1)
  extract(sys.argv[1], sys.argv[2])
