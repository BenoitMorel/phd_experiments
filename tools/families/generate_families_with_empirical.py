import sys
import os
import shutil
import ete3

def mymakedirs(path):
  try:
    os.makedirs(path)
  except:
    pass

def read_tree(tree):
  lines = open(tree).readlines()
  for line in lines:
    if (not line.startswith(">")):
      return ete3.Tree(line, format=1)
  return None

def generate_mapping_file(input_tree_file, mapping_file, treerecs_mapping_file):
  tree = read_tree(input_tree_file)
  leaves = tree.get_leaves()
  writer = open(mapping_file, "w")
  treerecs_writer = open(treerecs_mapping_file, "w")
  species_to_genes = {}
  for leaf in leaves:
    species = leaf.name.split("_")[0]
    gene = leaf.name
    if (not species in species_to_genes):
      species_to_genes[species] = []
    species_to_genes[species].append(gene)
  for species in species_to_genes:
    genes = species_to_genes[species]
    writer.write(species + ":")
    writer.write(";".join(genes))
    writer.write("\n")
    for gene in genes:
      treerecs_writer.write(species + " " + gene + "\n")

def generate_families_with_empirical(msas_dir, trees_dir, species_tree, dataset_dir):
  families_dir = os.path.join(dataset_dir, "families")
  all_alignments_dir = os.path.join(dataset_dir, "alignments")
  mymakedirs(families_dir)
  mymakedirs(all_alignments_dir)
  shutil.copy(species_tree, os.path.join(dataset_dir, "speciesTree.newick"))
  trees_dico = {}
  for tree in os.listdir(trees_dir):
    family_name = tree.split(".")[0]
    trees_dico[family_name] = os.path.join(trees_dir, tree)
  
  for msa in os.listdir(msas_dir):
    family_name = msa.split(".")[0]
    family_path = os.path.join(families_dir, family_name)
    mymakedirs(family_path)
    msa_source = os.path.join(msas_dir, msa)
    msa_dest = os.path.join(family_path, "alignment.msa")
    tree_source = trees_dico[family_name]
    tree_dest = os.path.join(family_path, "trueGeneTree.newick")
    mapping_file = os.path.join(family_path, "mapping.link")
    treerecs_mapping_file = os.path.join(family_path, "treerecs_mapping.link")
    shutil.copy(msa_source, msa_dest)
    shutil.copy(tree_source, tree_dest)
    shutil.copy(msa_source, os.path.join(all_alignments_dir, msa))
    generate_mapping_file(tree_source, mapping_file, treerecs_mapping_file)
    shutil.copy(species_tree, os.path.join(family_path, "speciesTree.newick"))



print("For now, this script should only work for cyano dataset")
if (len(sys.argv) != 5):
  print("Syntax: python generate_families_with_empirical.py msas_dir trees_dir species_tree dataset_diri")
  exit(1)

msas_dir = sys.argv[1]
trees_dir = sys.argv[2]
species_tree = sys.argv[3]
dataset_dir = sys.argv[4]
generate_families_with_empirical(msas_dir, trees_dir, species_tree, dataset_dir)


