import sys
import os
import shutil
import ete3
import fam
sys.path.insert(0, 'tools/msa_edition')
import msa_converter


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
      treerecs_writer.write(gene + " " + species + "\n")

def generate_families_with_empirical(msas_dir, trees_dir, species_tree, datadir, do_convert):
  fam.init_top_directories(datadir) 
  shutil.copy(species_tree, fam.get_species_tree(datadir))
  trees_dico = {}
  
  for tree in os.listdir(trees_dir):
    family_name = tree.split(".")[0]
    trees_dico[family_name] = os.path.join(trees_dir, tree)
  
  for msa in os.listdir(msas_dir):
    family = msa.split(".")[0]
    fam.init_family_directories(datadir, family)

    msa_source = os.path.join(msas_dir, msa)
    msa_dest = fam.get_alignment(datadir, family)
    tree_source = trees_dico[family]
    tree_dest = fam.get_true_tree(datadir, family)
    mapping_file = fam.get_mappings(datadir, family)
    treerecs_mapping_file = fam.get_treerecs_mappings(datadir, family)
    if (do_convert):
      msa_converter.msa_convert(msa_source, msa_dest, "phylip_relaxed", "fasta")
    else:
      shutil.copy(msa_source, msa_dest)
    shutil.copy(tree_source, tree_dest)
    generate_mapping_file(tree_source, mapping_file, treerecs_mapping_file)
  fam.postprocess_datadir(datadir)


if (len(sys.argv) != 6):
  print("Syntax: python generate_families_with_empirical.py msas_dir trees_dir species_tree datadiri do_convert_msas")
  print("set do_convert_msas to 1 for empirical and 0 for simulated")
  exit(1)

msas_dir = sys.argv[1]
trees_dir = sys.argv[2]
species_tree = sys.argv[3]
datadir = sys.argv[4]
do_convert = int(sys.argv[5]) != 0
generate_families_with_empirical(msas_dir, trees_dir, species_tree, datadir, do_convert)


