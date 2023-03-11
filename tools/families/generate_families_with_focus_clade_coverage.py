import sys
import os
import ete3
import math 
import shutil
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
sys.path.insert(0, 'tools/mappings')
import fam_data
import fam
import get_dico
sys.path.insert(0, 'scripts')
import experiments as exp
from read_tree import read_tree

def filter_ok(input_datadir, family, taxa, min_coverage):
  species_dict = get_dico.get_species_to_genes_family(input_datadir, family)
  count = 0
  for taxon in taxa:
    if (taxon in species_dict):
      count += 1
  return count >= min_coverage
  
def get_focus_taxa(species_tree, subtree_name):
  tree = read_tree(species_tree)
  subtree = None
  for node in tree.traverse():
    if (node.name == subtree_name):
      subtree = node
      break
  assert(subtree != None)
  return node.get_leaf_names()
  

def generate(input_datadir, species_tree, subtree, min_coverage):
  output_datadir = os.path.normpath(input_datadir)
  output_datadir += "_" + subtree
  output_datadir += "_cov" + str(min_coverage)
  if (os.path.exists(output_datadir)):
    print("Error, output datadir already exists " + output_datadir)
    sys.exit(1)
  taxa = get_focus_taxa(species_tree, subtree)
  print(taxa)
  fam.init_top_directories(output_datadir)   
  exp.relative_symlink(species_tree, fam.get_species_tree(output_datadir))
  families = fam.get_families_list(input_datadir)
  for family in families:
    if (filter_ok(input_datadir, family, taxa, min_coverage)):
      fam_data.duplicate_families_symlink(input_datadir, output_datadir, family)
  fam.postprocess_datadir(output_datadir)
  try:
    shutil.copyfile(fam.get_species_dict(input_datadir), fam.get_species_dict(output_datadir))
  except:
    pass
  output_families = fam.get_families_list(output_datadir)
  print("Input families: \t" + str(len(families)))
  print("Output families: \t" + str(len(output_families)))
  print("Done. Resulting datadir in " + output_datadir)
  return output_datadir



if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax python " + os.path.basename(__file__) + " datadir species_tree subtree min_coverage")
    sys.exit(1)

  input_datadir = sys.argv[1]
  species_tree = sys.argv[2]
  subtree = sys.argv[3]
  min_coverage = int(sys.argv[4])
  generate(input_datadir, species_tree, subtree, min_coverage)

