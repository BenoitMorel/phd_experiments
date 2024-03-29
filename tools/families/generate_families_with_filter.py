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

def filter_ok(datadir, family, min_species, max_species, min_sites):
  species_dict = get_dico.get_species_to_genes_family(datadir, family)
  if (len(species_dict) < min_species):
    return False
  if (len(species_dict) > max_species):
    return False
  if (min_sites > 0):
    msa_file = fam.get_alignment(datadir, family)
    seqs = ete3.SeqGroup(open(msa_file).read())
    first_seq = seqs.get_entries()[0][1]
    sites = len(first_seq)
    if (sites < min_sites):
      return False
  return True 

def get_output_dir(input_datadir, min_coverage_ratio, max_coverage_ratio, min_sites):
  output_datadir = os.path.normpath(input_datadir)
  output_datadir += "_mincov" + str(min_coverage_ratio) 
  output_datadir += "_maxcov" + str(max_coverage_ratio) 
  if (min_sites > 0):
    output_datadir += "_minsites" + str(min_sites)
  return output_datadir

def generate(input_datadir, min_coverage_ratio, max_coverage_ratio, min_sites):
  output_datadir = get_output_dir(input_datadir, min_coverage_ratio, max_coverage_ratio, min_sites)
  if (os.path.exists(output_datadir)):
    print("Error, output datadir already exists " + output_datadir)
    sys.exit(1)
  species_tree = fam.get_species_tree(input_datadir)
  species = len(ete3.Tree(species_tree, 1).get_leaves())
  min_species = int(float(species) * min_coverage_ratio)
  max_species = int(math.ceil(float(species) * max_coverage_ratio))
  fam.init_top_directories(output_datadir)   
  exp.relative_symlink(species_tree, fam.get_species_tree(output_datadir))
  families = fam.get_families_list(input_datadir)
  for family in families:
    if (filter_ok(input_datadir, family, min_species, max_species, min_sites)):
      fam_data.duplicate_families_symlink(input_datadir, output_datadir, family)
  fam.postprocess_datadir(output_datadir)
  try:
    shutil.copyfile(fam.get_species_dict(input_datadir), fam.get_species_dict(output_datadir))
  except:
    pass
  output_families = fam.get_families_list(output_datadir)
  print("Min families " + str(min_species))
  print("Input families: \t" + str(len(families)))
  print("Output families: \t" + str(len(output_families)))
  print("Done. Resulting datadir in " + output_datadir)
  return output_datadir



if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax python " + os.path.basename(__file__) + " datadir min_coverage_ratio max_coverage_ratio min_sites")
    sys.exit(1)

  input_datadir = sys.argv[1]
  min_coverage_ratio = float(sys.argv[2])
  max_coverage_ratio = float(sys.argv[3])
  min_sites = int(sys.argv[4])
  generate(input_datadir, min_coverage_ratio, max_coverage_ratio, min_sites)
