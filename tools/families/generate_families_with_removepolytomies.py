import sys
import os
import ete3
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/mappings')
import fam_data
import fam
import get_dico
sys.path.insert(0, 'scripts')
import experiments as exp

def filter_ok(datadir, family):
  msa_file = fam.get_alignment(datadir, family)
  seqs = ete3.SeqGroup(open(msa_file).read())
  counter = {}
  for entry in seqs.get_entries():
    seq = entry[1]
    if (seq in counter):
      counter[seq] += 1
      if (counter[seq] > 2):
        return False
    else:
      counter[seq] = 1
  return True

def generate(input_datadir):
  output_datadir = os.path.normpath(input_datadir)
  output_datadir += "_nopoly"
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
  output_families = fam.get_families_list(output_datadir)
  print("Input families: \t" + str(len(families)))
  print("Output families: \t" + str(len(output_families)))
  print("Done. Resulting datadir in " + output_datadir)




if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " datadir")
    sys.exit(1)

  input_datadir = sys.argv[1]
  generate(input_datadir)

