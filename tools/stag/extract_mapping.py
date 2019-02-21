import os
import sys
from ete3 import Tree


def extract_from_species_tree(species_tree_path, output):
  tree = Tree(open(species_tree_path).read(), format=1) 
  with open(output, "w") as writer:
    for leaf in tree.get_leaves():
      writer.write(leaf.name + "_* " + leaf.name + "\n")

def extract_from_dataset(dataset, output):
  dico = {}
  families_path = os.path.join(dataset, "families")
  species_tree = Tree(open(os.path.join(dataset, "speciesTree.newick")).read(), format=1)
  for leaf in species_tree.get_leaves():
    dico[leaf.name] = []
  for family in os.listdir(families_path):
    mapping_file = os.path.join(families_path, family, "mapping.link")
    for line in open(mapping_file).readlines():
      sp1 = line.split(":")
      species = sp1[0]
      genes = sp1[1][:-1].split(";")
      dico[species].extend(genes)
  
  with open(output, "w") as writer:
    for species in dico:
      for gene in dico[species]:
        writer.write(gene + " " + species + "\n")
      #writer.write(species + " " + " ".join(dico[species]) + "\n")


if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax: python extract_mapping.py [species_tree, dataset] output_mapping")
    exit(1)

  if (os.path.isfile(sys.argv[1])):
    extract_from_species_tree(sys.argv[1], sys.argv[2])
  else:
    extract_from_dataset(sys.argv[1], sys.argv[2])



