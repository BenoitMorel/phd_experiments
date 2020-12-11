import sys
import os
import ete3

def get_taxa_number(msa_file):
  seqs = ete3.SeqGroup(open(msa_file).read()) #, format="phylip_relaxed")
  return len(seqs.get_entries())

if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " msa")
    sys.exit(1)
  print("Taxa number: " + str(get_taxa_number(sys.argv[1])))
