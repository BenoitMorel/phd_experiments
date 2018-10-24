import os
import sys
from ete3 import Tree


if (len(sys.argv) != 2):
  print("syntax: families_dir")
  exit(1)


families_dir = sys.argv[1]


def generate_link_file(gene_tree_path, output_link_path):
  print(gene_tree_path)
  if (len(open(gene_tree_path).readlines()) == 0):
    return
  gene_tree = Tree(gene_tree_path, format=1)

  leaves = gene_tree.get_leaf_names()

  dico = {}
  for leaf in leaves:
    split = leaf.split("@")
    species = split[0]
    if (not (species in dico)):
      dico[species] = []
    dico[species].append(leaf)

  with open(output_link_path, "w") as writer:
    for key, value in dico.items():
      writer.write(key)
      writer.write(":")
      writer.write(";".join(value))
      writer.write("\n")

for family in os.listdir(families_dir):
  family_dir = os.path.join(families_dir, family)
  gene_tree_path = os.path.join(family_dir, "raxmlGeneTree.newick")
  output_link_path = os.path.join(family_dir, "phyldogMapping.link")
  generate_link_file(gene_tree_path, output_link_path)
