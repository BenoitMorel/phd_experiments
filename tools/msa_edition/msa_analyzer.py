import sys
import os
import ete3

def get_taxa_number(msa_file):
  seqs = ete3.SeqGroup(open(msa_file).read()) #, format="phylip_relaxed")
  return len(seqs.get_entries())


