import sys
import os
import shutil
import functools
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
from ete3 import Tree

# download data from https://github.com/chaoszhang/A-pro_data


def extract(gene_trees, species_tree, datadir, cut_underscore = False):
  fam.init_top_directories(datadir)
  index = 0
  if (species_tree != None):
    shutil.copyfile(species_tree, fam.get_species_tree(datadir))
    print("\tCopied species tree into " + fam.get_species_tree(datadir))
  for line in open(gene_trees).readlines():
    if (line[0] == "#"):
      continue
    family = "family_" + str(index)
    tree_str = line.replace("[&U]", "")
    tree = Tree(tree_str)
    leaves = tree.get_leaves()
    if (len(leaves) < 4):
      continue
    leaves_dict = {}
    species_to_genes = {}
    for leaf in leaves:
      species = leaf.name
      if (cut_underscore):
        species = species.split("_")[0]
      if (not species in species_to_genes):
        species_to_genes[species] = []
      genes = species_to_genes[species]
      leaf.name = species + "_" + str(len(genes))
      genes.append(leaf.name)
    fam.init_family_directories(datadir, family)
    mapping_file = fam.get_mappings(datadir, family)
    fam.write_phyldog_mapping(species_to_genes, mapping_file)
    with open(fam.get_true_tree(datadir, family), "w") as writer:
      writer.write(tree.write())
    index += 1
  print("\tPostprocessing data...")
  fam.postprocess_datadir(datadir)

def extract_fungi(inputdir, datadir):
  gene_trees_path = os.path.join(inputdir, "pep.ml.renamed.trees")
  species_tree = os.path.join(inputdir, "fungi-ref-nonewline.tre")
  print("Extracting fungi astral-pro data...")
  extract(gene_trees_path, species_tree, datadir)

def extract_plants(inputdir, datadir):
  gene_trees_path = os.path.join(inputdir, "1kp-c12-genetrees.tre")
  astral_species_tree = os.path.join(inputdir, "single-copy-astral-localPP-recomputed-originalformula.tre")
  print("Extracting 1kplant astral-pro data...")
  extract(gene_trees_path, astral_species_tree, datadir)

def extract_tom(datadir):
  inputdir = os.path.join(exp.benoit_datasets_root, "raw_data", "life92")
  gene_trees_path = os.path.join(inputdir, "gene_trees.txt")
  species_tree = None
  extract(gene_trees_path, species_tree, datadir)

def extract_bigcyano36(datadir):
  inputdir = os.path.join(exp.benoit_datasets_root, "raw_data", "bigcyano")
  gene_trees_path = os.path.join(inputdir, "gene_trees.txt")
  species_tree = None
  extract(gene_trees_path, species_tree, datadir, True)


def extract_aprodata(rawdatadir):
  families_dir = os.path.join(exp.benoit_datasets_root, "families")
  fungi_input_dir = os.path.join(rawdatadir, "fungi")
  fungi_output_dir = os.path.join(families_dir, "apro_fungi")
  extract_fungi(fungi_input_dir, fungi_output_dir)
  plants_input_dir = os.path.join(rawdatadir, "1kp")
  plants_output_dir = os.path.join(families_dir, "apro_plants")
  #extract_plants(plants_input_dir, plants_output_dir)



if (__name__ == "__main__"): 
  if (len(sys.argv) < 2): 
    print("Syntax: python " + os.path.basename(__file__) + " astralprodata_repository_path")
    exit(1)
  #extract_aprodata(sys.argv[1])
  extract_tom(sys.argv[1])
