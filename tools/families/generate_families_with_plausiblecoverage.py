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

def get_families_list(input_datadir):
  input_families = fam.get_families_list(input_datadir)
  output_families = []
  for family in input_families:
    coverage = []
    species_to_genes = get_dico.get_species_to_genes_family(input_datadir, family)
    found_dup = False
    total_genes = 0
    for species in species_to_genes:
      genes = len(species_to_genes[species])
      coverage.append(genes)
      total_genes += genes
    av = float(sum(coverage)) / float(len(coverage))
    if (family == "family_1310"):
      print(coverage)
    if (max(coverage) < 2 * av):
      output_families.append(family)
    else:
      print("Remove family " + family + " with " + str(total_genes) + " genes")
  print("Input families: " + str(len(input_families)))
  print("Output families: " + str(len(output_families)))
  return output_families


def generate(input_datadir):
  output_datadir = input_datadir + "_plausible_coverage"
  if (os.path.exists(output_datadir)):
    print("Directory " + output_datadir + " already exists. Skipping.")
    return
  
  fam.init_top_directories(output_datadir)   
  species_tree = fam.get_species_tree(input_datadir)
  exp.relative_symlink(species_tree, fam.get_species_tree(output_datadir))
  families = get_families_list(input_datadir)
  for family in families:
    fam_data.duplicate_families_symlink(input_datadir, output_datadir, family)
  fam.postprocess_datadir(output_datadir)
  output_families = fam.get_families_list(output_datadir)
  print("Result datadir in " + output_datadir)




if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " datadir")
    sys.exit(1)

  input_datadir = sys.argv[1]
  generate(input_datadir)

