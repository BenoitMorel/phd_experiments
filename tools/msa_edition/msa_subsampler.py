import os
import sys
import ete3
from random import shuffle



def subsample(input_msa, output_msa, msa_format, ratio):
  seqs = ete3.SeqGroup(input_msa, format=msa_format)
  sites = 0
  sites = len(seqs.get_entries()[0][1])
  indices = range(sites)
  shuffle(indices)
  new_sites = int(float(sites) * ratio)
  indices = indices[0: new_sites]
  for seq in  seqs.get_entries():
    new_seq = ""
    for i in indices:
      new_seq += seq[1][i]
    seqs.set_seq(seq[0], new_seq)
  seqs.write(msa_format, output_msa)


"""
  Remove all sequences that are not in to_keep
"""
def prune_sequences(input_msa, output_msa, msa_format, to_keep):
  seqs = ete3.SeqGroup(open(input_msa).read(), format=msa_format)
  new_msa = ete3.SeqGroup()
  for entry in seqs.get_entries():
    if (entry[0] in to_keep):
      new_msa.set_seq(entry[0], entry[1])
  new_msa.write(msa_format, output_msa)

if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax: python msa_subsampler input_msa output_msa format ratio")
    print("Format are: fasta, phylip, iphylip, phylip_relaxed, iphylip_relaxed")
    exit(1)
    
  input_msa = sys.argv[1]
  output_msa = sys.argv[2]
  msa_format = sys.argv[3]
  ratio = float(sys.argv[4])
  subsample(input_msa, output_msa, msa_format, ratio)
