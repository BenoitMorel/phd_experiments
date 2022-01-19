import sys
import os
from ete3 import Tree
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/mappings')
import experiments as exp
import fam
import get_dico

def get_species_dict(datadir, species_tree_path, dictionary_path):
  d = {}
  tree = Tree(species_tree_path, format = 1)
  for leaf in tree.get_leaves():
    d[leaf.name] = leaf.name
  if (dictionary_path == None):
    dictionary_path = fam.get_species_dict(datadir)
  print(dictionary_path)
  if (not os.path.isfile(dictionary_path)):
    print("No datadir dict found")
    return d
  for line in open(dictionary_path).readlines():
    sp = line.replace("\n", "").split(":")
    d[sp[0]] = sp[1]
  return  d

def get_unique(name, existing_names):
  i = 0
  while (True):
    unique_name = name + "_" + str(i)
    if (not unique_name in existing_names):
      existing_names.add(unique_name)
      return unique_name
    i += 1


def translate(gene_tree_path, dictionary_path):
  tree = Tree(gene_tree_path, format = 1)
  sp = os.path.normpath(gene_tree_path).split(os.sep)
  family = sp[-3]
  datadir = gene_tree_path.split(os.path.join(sp[-4], sp[-3]))[0]
  species_tree_path = fam.get_species_tree(datadir) 
  gene_to_species = get_dico.get_gene_to_species(datadir, family)
  existing_names = set()
  species_dict = get_species_dict(datadir, gene_tree_path, dictionary_path)
  for leaf in tree.get_leaves():
    species = gene_to_species[leaf.name]
    trans_species = species_dict[species]
    leaf.name = get_unique(trans_species, existing_names)
  return tree.write(format = 1)

def dump_into(gene_tree_path, output_path):
  with open(output_path, "w") as writer:
    writer.write(translate(gene_tree_path))


if (__name__ == "__main__"): 
  if (len(sys.argv) < 2): 
    print("Syntax: python " + os.path.basename(__file__) + " gene_tree [dictionary]")
    exit(1)

  gene_tree_path = sys.argv[1]
  dictionary_path = None
  if (len(sys.argv) > 2):
    dictionary_path = sys.argv[2]
  print(translate(gene_tree_path, dictionary_path))


