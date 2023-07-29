import sys
import os
import ete3
import math 
import shutil
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/mappings')
import fam_data
import fam
import get_dico
sys.path.insert(0, 'scripts')
import experiments as exp

def filter_ok(datadir, family):
  species_dict = get_dico.get_species_to_genes_family(datadir, family)
  for species in species_dict:
    if (len(species_dict[species]) > 1):
      return True
  return False


def generate(input_datadir):
  output_datadir = input_datadir + "_multicopy"
  if (os.path.exists(output_datadir)):
    print("Error, output datadir already exists " + output_datadir)
    sys.exit(1)
  species_tree = fam.get_species_tree(input_datadir)
  fam.init_top_directories(output_datadir)   
  exp.relative_symlink(species_tree, fam.get_species_tree(output_datadir))
  families = fam.get_families_list(input_datadir)
  for family in families:
    if (filter_ok(input_datadir, family)):
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
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " datadir")
    sys.exit(1)

  input_datadir = sys.argv[1]
  generate(input_datadir)

