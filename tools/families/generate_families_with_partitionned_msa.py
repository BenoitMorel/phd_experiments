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


is_dna = False

def generate(msa_file, partition_file, species_tree_file, datadir):
  alignments = splitter.get_sub_alignments(msa_file, partition_file)
  fam.init_top_directories(datadir)
  
  total = 0
  mini = 11198
  maxi = 0
  species = set()
  for family in alignments:
    alignment = remove_empty_sequences.get_cleaned_msa(alignments[family], is_dna)
    
    genes = {}
    for seq in alignment.get_entries():
      genes[seq[0]] = [seq[0]]
      species.add(seq[0])
    if (len(genes) < 4):
      continue
    total += len(genes)
    mini = min(mini, len(genes))
    maxi = max(maxi, len(genes))
    fam.init_family_directories(datadir, family)
    mapping_file = fam.get_mappings(datadir, family)
    fam.write_phyldog_mapping(genes, mapping_file)
    output_alignment = fam.get_alignment(datadir, family)
    alignment.write(format = "fasta", outfile = output_alignment)
  fam.postprocess_datadir(datadir)
  
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
  if (len(sys.argv) < 5): 
    print("Syntax: python " + os.path.basename(__file__) + " msa partition species_tree datadir")
    exit(1)
  msa_file = sys.argv[1]
  partition_file = sys.argv[2]
  species_tree_file = sys.argv[3]
  datadir = sys.argv[4]
  generate(msa_file, partition_file, species_tree_file, datadir)

