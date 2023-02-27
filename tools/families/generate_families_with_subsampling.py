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

def check_families_coverage(input_datadir, families):
  species_tree = fam.get_species_tree(input_datadir)
  species_set = set(ete3.Tree(species_tree, 1).get_leaf_names())
  new_species_set = set()
  for family in families:
    for species in get_dico.get_species_to_genes_family(input_datadir, family):
      new_species_set.add(species)
  return new_species_set == species_set

def get_families_list(input_datadir, sampling_ratio):
  input_families = fam.get_families_list(input_datadir)
  output_families_number = int(float(len(input_families)) * sampling_ratio)
  for trial in range(0, 10):
    output_families = random.sample(input_families, output_families_number)
    if (check_families_coverage(input_datadir, output_families)):
      return output_families
  require_full_coverage = False
  print("Failed to find a sample of families covering all species")
  if (require_full_coverage):
    sys.exit(1)
  else:
    return random.sample(input_families, output_families_number)



def generate_replicate(input_datadir, sampling_ratio, replicate):
  random.seed(replicate + 42)
  input_datadir = os.path.normpath(input_datadir)
  output_datadir = input_datadir + "_subsample" + str(sampling_ratio)
  output_datadir += "_rep" + str(replicate) 
  if (os.path.exists(output_datadir)):
    print("Directory " + output_datadir + " already exists. Skipping.")
    return
  
  fam.init_top_directories(output_datadir)   
  species_tree = fam.get_species_tree(input_datadir)
  exp.relative_symlink(species_tree, fam.get_species_tree(output_datadir))
  families = get_families_list(input_datadir, sampling_ratio)
  for family in families:
    fam_data.duplicate_families_symlink(input_datadir, output_datadir, family)
  fam.postprocess_datadir(output_datadir)
  output_families = fam.get_families_list(output_datadir)
  print("Result datadir in " + output_datadir)


def generate(input_datadir, sampling_ratio, replicates):
  for replicate in replicates:
    generate_replicate(input_datadir, sampling_ratio, replicate)




if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python " + os.path.basename(__file__) + " datadir sampling_ratio replicates")
    sys.exit(1)

  input_datadir = sys.argv[1]
  sampling_ratio = float(sys.argv[2])
  replicates = int(sys.argv[3])
  generate(input_datadir, sampling_ratio, range(0, replicates))

