import sys
import os
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
sys.path.insert(0, 'tools/msa_edition')
import experiments as exp
import fam
import ete3
import read_msa
import split_partitionned_alignment as splitter 
import remove_empty_sequences
import create_random_tree



def generate(msa_file, partition_file, is_dna, species_tree_file, datadir):
  alignments = splitter.get_sub_alignments(msa_file, partition_file)
  fam.init_top_directories(datadir)
  
  total = 0
  mini = 999999999
  maxi = 0
  species = set()
  for family in alignments:
    alignment = remove_empty_sequences.get_cleaned_msa(alignments[family], is_dna)
    
    genes = {}
    for seq in alignment.get_entries():
      genes[seq[0]] = [seq[0]]
    print(len(genes))
    if (len(genes) < 4):
      continue
    for seq in alignment.get_entries():
      species.add(seq[0])
    total += len(genes)
    mini = min(mini, len(genes))
    maxi = max(maxi, len(genes))
    fam.init_family_directories(datadir, family)
    mapping_file = fam.get_mappings(datadir, family)
    fam.write_phyldog_mapping(genes, mapping_file)
    output_alignment = fam.get_alignment(datadir, family)
    alignment.write(format = "fasta", outfile = output_alignment)
    print(family)
  fam.postprocess_datadir(datadir)
  print("Number of species: " + str(len(species))) 
  true_species_tree = fam.get_species_tree(datadir)
  try:
    shutil.copy(species_tree_file, true_species_tree)
  except:
    print("Incorrect input species tree, creating a random tree")
    tree = create_random_tree.create_random_tree_from_species(species)
    tree.write(format= 1, outfile = true_species_tree)
  print(mini)
  print(maxi)
  print(float(total) / float(len(alignments)))


if (__name__ == "__main__"): 
  if (len(sys.argv) < 6): 
    print("Syntax: python " + os.path.basename(__file__) + " msa partition is_dna species_tree datadir")
    exit(1)
  msa_file = sys.argv[1]
  partition_file = sys.argv[2]
  is_dna = (0 != int(sys.argv[3]))
  species_tree_file = sys.argv[4]
  datadir = sys.argv[5]
  generate(msa_file, partition_file, is_dna, species_tree_file, datadir)

