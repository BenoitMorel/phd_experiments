import sys
import os
import shutil
import functools
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
import ete3

def get_alignments(rawdir):
  alignments_dir = os.path.join(rawdir, "alns")
  per_family_alignments = {}
  for alignments_subdir in os.listdir(alignments_dir):
    if (alignments_subdir.startswith(".")):
      continue
    alignments_subdir = os.path.join(alignments_dir, alignments_subdir)
    for alignment in os.listdir(alignments_subdir):
      family = alignment.replace(".fasta", "").replace(".", "_")
      per_family_alignments[family] = os.path.join(alignments_subdir, alignment)
  return per_family_alignments

def get_trees(rawdir):
  pargenes_dir = os.path.join(rawdir, "trees")
  per_family_trees = {}
  for pargenes_subdir in os.listdir(pargenes_dir):
    if (pargenes_subdir.startswith(".")):
      continue
    pargenes_subdir = os.path.join(pargenes_dir, pargenes_subdir)
    results = os.path.join(pargenes_subdir, "mlsearch_run", "results")
    for result in os.listdir(results):
      if (result.startswith(".")):
        continue
      family = result.replace("_fasta", "").replace(".", "_")
      per_family_trees[family] = os.path.join(results, result, result + ".raxml.bestTree")
  return per_family_trees 

def create_mapping_from_gene_tree(input_tree, output_mapping):
  leaves = ete3.Tree(input_tree, format=1).get_leaf_names()
  with open(output_mapping, "w") as writer:
    for leaf in leaves:
      species = leaf.split("_")[1]
      writer.write(leaf + ":" + species + "\n")


    

def get_families_from_trees(per_family_trees):
  families = []
  for family in per_family_trees:
    families.append(family)
  return families

def generate(rawdir, species_tree, datadir):
  
  # init directories
  print("Starts generation")
  fam.init_top_directories(datadir)
  
  # species tree
  true_species_tree = fam.get_species_tree(datadir)
  shutil.copy(species_tree, true_species_tree)
  fam.init_top_directories(datadir)
  
  # families 
  print("Init families")
  per_family_alignments = get_alignments(rawdir)
  per_family_trees = get_trees(rawdir)
  families = get_families_from_trees(per_family_trees)
  fam.init_families_directories(datadir, families)

  # fill families
  print("Fillf amilies")
  for family in families:
    # alignments
    src_alignment = per_family_alignments[family]
    dest_alignment = fam.get_alignment(datadir, family)
    shutil.copyfile(src_alignment, dest_alignment)

    # trees
    src_tree = per_family_trees[family]
    dest_tree = fam.get_true_tree(datadir, family)
    shutil.copyfile(src_tree, dest_tree)

    # mappings
    create_mapping_from_gene_tree(dest_tree, fam.get_mappings(datadir, family))

  print("post process")
  fam.postprocess_datadir(datadir)

if (__name__ == "__main__"): 
  if (len(sys.argv) < 4): 
    print("Syntax: python " + os.path.basename(__file__) + " raw_dir species_tree datadir")
    exit(1)
  rawdir = sys.argv[1]
  species_tree = sys.argv[2]
  datadir = sys.argv[3]
  generate(rawdir, species_tree, datadir)


