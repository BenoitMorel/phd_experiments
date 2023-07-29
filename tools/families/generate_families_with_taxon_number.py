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

def filter_ok(datadir, family, min_taxa, max_taxa):
  genes_dict = get_dico.get_gene_to_species(datadir, family)
  taxa = len(genes_dict)
  return taxa >= min_taxa and taxa <= max_taxa

def get_output_dir(input_datadir, min_taxa, max_taxa):
  output_datadir = os.path.normpath(input_datadir)
  if (min_taxa >= 0):
    output_datadir += "_mintaxa" + str(min_taxa) 
  if (max_taxa >= 0):
    output_datadir += "_maxtaxa" + str(max_taxa) 
  return output_datadir

def generate(input_datadir, min_taxa, max_taxa):
  output_datadir = get_output_dir(input_datadir, min_taxa, max_taxa)
  if (os.path.exists(output_datadir)):
    print("Error, output datadir already exists " + output_datadir)
    sys.exit(1)
  species_tree = fam.get_species_tree(input_datadir)
  fam.init_top_directories(output_datadir)   
  exp.relative_symlink(species_tree, fam.get_species_tree(output_datadir))
  families = fam.get_families_list(input_datadir)
  for family in families:
    if (filter_ok(input_datadir, family, min_taxa, max_taxa)):
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
  if (len(sys.argv) != 4):
    print("Syntax python " + os.path.basename(__file__) + " datadir min_taxa max_taxa")
    sys.exit(1)

  input_datadir = sys.argv[1]
  min_taxa = int(sys.argv[2])
  max_taxa = int(sys.argv[3])
  generate(input_datadir, min_taxa, max_taxa)

