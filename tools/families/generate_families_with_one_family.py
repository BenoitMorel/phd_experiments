import sys
import os
import shutil
sys.path.insert(0, 'tools/families')
import fam_data
import fam
sys.path.insert(0, 'scripts')
import experiments as exp



def generate(input_datadir, family):
  input_datadir = os.path.normpath(input_datadir)
  output_datadir = input_datadir + "_" + family
  if (os.path.exists(output_datadir)):
    print("Error, output datadir already exists " + output_datadir)
    sys.exit(1)
  species_tree = fam.get_species_tree(input_datadir)
  fam.init_top_directories(output_datadir)   
  shutil.copyfile(species_tree, fam.get_species_tree(output_datadir))
  fam_data.duplicate_families(input_datadir, output_datadir, family)
  fam.postprocess_datadir(output_datadir)
  try:
    shutil.copyfile(fam.get_species_dict(input_datadir), fam.get_species_dict(output_datadir))
  except:
    pass
  return output_datadir


if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python " + os.path.basename(__file__) + " datadir family")
    sys.exit(1)

  input_datadir = sys.argv[1]
  family = sys.argv[2]
  generate(input_datadir, family)


