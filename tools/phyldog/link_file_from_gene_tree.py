import os
import sys
from ete3 import Tree


if (len(sys.argv) != 3):
  print("syntax: input_gene_tree output_link_file")
  exit(1)

gene_tree_path = sys.argv[1]
output_link_path = sys.argv[2]

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


