import sys
import os
import ete3
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/mappings')
import fam_data
import fam
import get_dico
import shutil
sys.path.insert(0, 'scripts')
import experiments as exp
from ete3 import SeqGroup



def filter_and_copy(input_datadir, output_datadir, family, is_dna):
  
  seqs = SeqGroup(fam.get_alignment(input_datadir, family))
  seq_to_labels = {}
  label_to_seq = {}
  for entry in seqs.iter_entries():
    seq = ""
    # for our purpose, unkown and deletion are the same
    if (is_dna):
      seq = entry[1].replace("N", "-")
    else:
      seq = entry[1].replace("X", "-")
    if (seq in seq_to_labels):
      seq_to_labels[seq].append(entry[0])
    else:
      seq_to_labels[seq] = [entry[0]]
    label_to_seq[entry[0]] = seq
  # will contain gene labels that have at most
  # one other gene with the exact same sequence
  labels_to_keep = set()
  for duplicates in seq_to_labels.values():
    if (len(duplicates) <= 2):
      for label in duplicates:
        labels_to_keep.add(label)


  gene_to_species = get_dico.get_gene_to_species(input_datadir, family)
  new_species_to_genes = {}
  for label in labels_to_keep:
    species = gene_to_species[label]
    if (species in new_species_to_genes):
      new_species_to_genes[species].append(label)
    else:
      new_species_to_genes[species] = [label]
  MIN_COVERED_SPECIES = 3
  # remove families that are not covered by enough families 
  if (len(new_species_to_genes) < MIN_COVERED_SPECIES):
    return
  # now let's copy!
  fam.init_family_directories(output_datadir, family)
  output_alignment = fam.get_alignment(output_datadir, family)
  new_seqs = SeqGroup()
  for label in labels_to_keep:
    new_seqs.set_seq(label, label_to_seq[label])
  new_seqs.write(outfile = output_alignment)
  get_dico.export_species_to_genes(new_species_to_genes, output_datadir, family)


def generate(input_datadir, is_dna):
  output_datadir = input_datadir + "_nodup"
  if (os.path.exists(output_datadir)):
    print("Error, output datadir already exists " + output_datadir)
    sys.exit(1)
  fam.init_top_directories(output_datadir)   
  species_tree = fam.get_species_tree(input_datadir)
  shutil.copy(species_tree, fam.get_species_tree(output_datadir))
  families = fam.get_families_list(input_datadir)
  for family in families:
    print("Treating family " + family)
    filter_and_copy(input_datadir, output_datadir, family, is_dna)
  fam.postprocess_datadir(output_datadir)
  output_families = fam.get_families_list(output_datadir)
  print("Input families: \t" + str(len(families)))
  print("Output families: \t" + str(len(output_families)))
  print("Done. Resulting datadir in " + output_datadir)
  return output_datadir

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python " + os.path.basename(__file__) + " datadir is_dna")
    print("1 for dna, 0 for aa")
    sys.exit(1)

  input_datadir = sys.argv[1]
  is_dna = int(sys.argv[2]) == 1
  generate(input_datadir, is_dna)

