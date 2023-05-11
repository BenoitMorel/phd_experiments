import sys
import os
import random
import ete3
import fam
sys.path.insert(0, 'scripts')
import generate_families_with_prunespecies
import experiments as exp



def generate_replicate(input_datadir, sampling_ratio, replicate):
  random.seed(replicate + 42)
  input_datadir = os.path.normpath(input_datadir)
  output_datadir = input_datadir + "_subtax" + str(sampling_ratio)
  output_datadir += "_rep" + str(replicate) 
  if (os.path.exists(output_datadir)):
    print("Directory " + output_datadir + " already exists. Skipping.")
    return
  species_tree = fam.get_species_tree(input_datadir) 
  leaves = ete3.Tree(species_tree).get_leaf_names()
  number_to_remove = int(float(len(leaves)) * (1.0 - sampling_ratio))
  leaves_to_remove = random.sample(leaves, number_to_remove)
  generate_families_with_prunespecies.generate(input_datadir, output_datadir, "true", "true", False, leaves_to_remove)
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

