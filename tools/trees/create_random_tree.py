import ete3
import sys
import os
sys.path.insert(0, os.path.join("tools", "msa_edition"))
import read_msa
  
random_alignment_format = "fasta"

def create_random_tree(msa_file, output_file):
  global random_alignment_format
  msa = read_msa.read_msa(msa_file)

  tree = ete3.Tree()
  tree.populate(len(msa.get_entries()))

  leaves = tree.get_leaves()
  index = 0

  for entry in msa.get_entries():
    leaves[index].add_feature("name", entry[0])
    index += 1

  tree.write(outfile=output_file, format=1)

def create_random_tree_from_species(species):
  tree = ete3.Tree()
  tree.populate(len(species))
  leaves = tree.get_leaves()
  index = 0
  for s in species:
    leaves[index].name = s
    index += 1
  return tree

def create_random_tree_from_tree(tree_file, output_file = None):
  tree = ete3.Tree(tree_file)
  species = tree.get_leaf_names()
  random_tree = create_random_tree_from_species(species)
  if (None != output_file):
    random_tree.write(outfile = output_file)
  return random_tree

def create_random_tree_taxa_number(taxa_number):
  species = []
  for i in range(0, taxa_number):
    species.append("hCoV-19/Australia/NSW14/2020|EPI_ISL_" + str(i) + "|2020-03-03")
    #species.append("HCOV-19_CHINA_2020012" + str(i) + "_EPI_ISL_" + str(i) )
  return create_random_tree_from_species(species)


if (__name__ == "__main__"): 
  if (len(sys.argv) != 3):
    print("Syntax: python create_random_tree.py msa output_random_tree")
    sys.exit(1)

  msa_file = sys.argv[1]
  output_file = sys.argv[2]
  create_random_tree(msa_file, output_file)
