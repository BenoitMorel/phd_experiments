import sys
import os
import shutil
import functools
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
import create_random_tree
import ete3


def treat_alignment(input_file, output_file, output_mapping_file, authorized_species):
  #seq = ete3.SeqGroup(input_file, format="iphylip_relaxed")
  print(input_file)
  seq = ete3.SeqGroup(input_file, format="fasta")
  filtered_seq = ete3.SeqGroup()
  species_to_genes = {}
  for entry in seq.iter_entries():
    species = entry[0].split("_")[0]
    if (not species in authorized_species):
      continue
    if (not species in species_to_genes):
      species_to_genes[species] = []
    species_to_genes[species].append(entry[0])
    filtered_seq.set_seq(entry[0], entry[1])
  filtered_seq.write("fasta", output_file)
  with open(output_mapping_file, "w") as writer:
    for species, genes in species_to_genes.items():
      writer.write(species + ":" + ";".join(genes) + "\n")


"""
  In msas_dir, the msas are fasta files, and each sequence is prefixed
  with a species name followed by an underscore. 
  This function builds a families direcotry structure on which I can
  run my benchmarks
  It removes all the sequences that are mapped to a species that
  is not in the species tree
"""
def generate_from_msas(msas_dir, species_tree, datadir):
  # init directories
  fam.init_top_directories(datadir)
  
  # species tree
  true_species_tree = fam.get_species_tree(datadir)
  shutil.copy(species_tree, true_species_tree)
  authorized_species = ete3.Tree(true_species_tree, format=1).get_leaf_names()

  families = []
  for f in os.listdir(msas_dir):
    families.append(f.split(".")[0])
  fam.init_families_directories(datadir, families)
  for f in os.listdir(msas_dir):
    family = f.split(".")[0]
    src = os.path.join(msas_dir, f)
    dest = fam.get_alignment(datadir, family)
    mapping_dest = fam.get_mappings(datadir, family)
    treat_alignment(src, dest, mapping_dest, authorized_species)
  fam.postprocess_datadir(datadir)

if (__name__ == "__main__"): 
  if (len(sys.argv) < 4): 
    print("Syntax: python " + os.path.basename(__file__) + " msas_dir species_tree datadir")
    exit(1)
  msas_dir = sys.argv[1]
  species_tree = sys.argv[2]
  datadir = sys.argv[3]
  generate_from_msas(msas_dir, species_tree, datadir)

