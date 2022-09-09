import sys
import os
import ete3

def get_taxa_number(msa_file):
  seqs = ete3.SeqGroup(open(msa_file).read()) #, format="phylip_relaxed")
  return len(seqs.get_entries())

def print_info(msa_file, form):
  seqs = ete3.SeqGroup(open(msa_file).read(), format=form)
  taxa = len(seqs.get_entries())
  sites = len(seqs.get_entries()[0][1])
  print("Taxa number: " + str(taxa))
  print("Sites number: " + str(sites))

if (__name__ == "__main__"):
  if (len(sys.argv) < 2):
    print("Syntax python " + os.path.basename(__file__) + " msa [format]")
    sys.exit(1)
  form = "fasta"
  if (len(sys.argv) > 2):
    form = sys.argv[2]
  print_info(sys.argv[1], form)

  
