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


def generate(msa_file, partition_file, species_tree_file, datadir):
  alignments = splitter.get_sub_alignments(msa_file, partition_file)
  fam.init_top_directories(datadir)
  true_species_tree = fam.get_species_tree(datadir)
  total = 0
  mini = 98
  maxi = 0
  for family in alignments:
    alignment = remove_empty_sequences.get_cleaned_msa_dna(alignments[family])
    
    genes = {}
    for seq in alignment.get_entries():
      genes[seq[0]] = seq[0]
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

