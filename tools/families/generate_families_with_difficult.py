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


def generate(input_datadir, min_score, max_score):
  input_datadir = os.path.normpath(input_datadir)
  output_datadir = input_datadir
  if (min_score > 0):
    output_datadir = output_datadir + "_mindiff" + str(min_score)
  if (max_score < 1):
    output_datadir = output_datadir + "_maxdiff" + str(max_score)

  if (os.path.exists(output_datadir)):
    print("Directory " + output_datadir + " already exists. Skipping.")
    return
  
  fam.init_top_directories(output_datadir)   
  species_tree = fam.get_species_tree(input_datadir)
  exp.relative_symlink(species_tree, fam.get_species_tree(output_datadir))
  for family in fam.get_families_list(input_datadir):
    score = fam.get_pythia_score(input_datadir, family)
    if (score < max_score and score > min_score):
      fam_data.duplicate_families_symlink(input_datadir, output_datadir, family)
  fam.postprocess_datadir(output_datadir)
  output_families = fam.get_families_list(output_datadir)
  print("Result datadir in " + output_datadir)




if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python " + os.path.basename(__file__) + " datadir min_score max_score")
    sys.exit(1)

  input_datadir = sys.argv[1]
  min_score = float(sys.argv[2])
  max_score = float(sys.argv[3])
  generate(input_datadir, min_score, max_score)


