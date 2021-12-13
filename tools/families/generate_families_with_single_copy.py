import sys
import os
import random
import ete3
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/mappings')
import fam_data
import fam
import get_dico
sys.path.insert(0, 'scripts')
import experiments as exp

def is_single_copy(input_datadir, family):
  d = get_dico.get_species_to_genes_family(input_datadir, family)
  for species in d:
    if (len(d[species]) > 1):
      return False
  return True

def get_single_copy_families(input_datadir):
  input_families = fam.get_families_list(input_datadir)
  print("Input families: " + str(len(input_families)))
  output_families = []
  for family in input_families:
    if (is_single_copy(input_datadir, family)):
      output_families.append(family)
  print("Single-copy families: " + str(len(output_families)))
  return output_families


def generate(input_datadir, output_datadir):
  if (os.path.exists(output_datadir)):
    print("Directory " + output_datadir + " already exists. Skipping.")
    return
  
  fam.init_top_directories(output_datadir)   
  species_tree = fam.get_species_tree(input_datadir)
  exp.relative_symlink(species_tree, fam.get_species_tree(output_datadir))
  families = get_single_copy_families(input_datadir)
  for family in families:
    fam_data.duplicate_families_symlink(input_datadir, output_datadir, family)
  fam.postprocess_datadir(output_datadir)
  output_families = fam.get_families_list(output_datadir)
  print("Result datadir in " + output_datadir)






if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " inputdir")
    sys.exit(1)

  input_datadir = sys.argv[1]
  output_datadir = os.path.normpath(input_datadir) + "_single"
  generate(input_datadir, output_datadir)




