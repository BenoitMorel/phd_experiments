import sys
import os
from Bio import SeqIO

def convert(nexus, fasta):
  records = SeqIO.parse(nexus, "nexus")
  count = SeqIO.write(records, fasta, "fasta")
  print("Converted %i records" % count)

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python " + os.path.basename(__file__) + " datadir gene_trees subst_model")
    sys.exit(1)
  nexus = sys.argv[1]
  fasta = sys.argv[2]
  convert(nexus, fasta)



