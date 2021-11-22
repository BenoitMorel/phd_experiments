import sys
import os
import shutil
import subprocess
import functools
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
import create_random_tree
from ete3 import Tree
from ete3 import SeqGroup
import re

def read_seqs(alignment_file):
  msa = None
  formats = ["fasta", "phylip_relaxed", "iphylip_relaxed", "phylip_interleaved"]
  for f in formats:
    try:
      msa = SeqGroup(alignment_file, f)
      return msa
    except:
      pass
  return None

def fullmatch(regex, string, flags=0):
  m = re.match("(?:" + regex + r")\Z", string, flags=flags)
  if (m):
    return True
  else:
    return False

def only_gaps(s):
  return fullmatch("[-?]*", s)

def prune_seqs(seqs):
  new_seqs = SeqGroup()
  for seq in  seqs.get_entries():
    if (not only_gaps(seq[1])):
      new_seqs.set_seq(seq[0], seq[1])
  return new_seqs

def generate(index, is_dna):
  input_dir = exp.stamatak_tests_dir 
  output_dir = os.path.join(exp.families_datasets_root, "stam")
  if (is_dna):
    input_dir = os.path.join(input_dir, "DNA-Data", str(index))
    output_dir = output_dir + "_DNA_" + str(index)
  else:
    input_dir = os.path.join(input_dir, "Protein-Data", str(index))
    output_dir = output_dir + "_AA_" + str(index)
  print("Input: " + input_dir)
  print("Output: " + output_dir)

  fam.init_top_directories(output_dir)
  partition_dir = os.path.join(input_dir, "split-partitions")
  species = set()
  for ali in os.listdir(partition_dir):
    family = ali.split(".")[0]
    print("Treating family " + family)
    fam.init_family_directories(output_dir, family)
    seqs = read_seqs(os.path.join(partition_dir, ali))
    print(os.path.join(partition_dir, ali))
    output_alignment = fam.get_alignment(output_dir, family)
    pruned_seqs = prune_seqs(seqs)
    pruned_seqs.write("fasta", output_alignment) 
    species_to_genes = {}
    for seq in pruned_seqs.get_entries():
      species.add(seq[0])
      species_to_genes[seq[0]] = [seq[0]]
    mapping_file = fam.get_mappings(output_dir, family)
    fam.write_phyldog_mapping(species_to_genes, mapping_file)

  random_species_tree = create_random_tree.create_random_tree_from_species(species)
  random_species_tree.write(outfile = fam.get_species_tree(output_dir))
  fam.postprocess_datadir(output_dir)
   

if (__name__ == "__main__"): 
  if (len(sys.argv) != 3): 
    print("Syntax: python " + os.path.basename(__file__) + " ID isDNA[1 for DNA and 0 for AA")
    exit(1)
  index = sys.argv[1]
  is_dna = sys.argv[2]
  assert(is_dna == "1" or is_dna == "0")
  is_dna = (is_dna == "1")
  generate(index, is_dna)


