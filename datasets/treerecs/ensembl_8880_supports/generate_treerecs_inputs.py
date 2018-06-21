import os 
import sys


fasta_dir = "/hits/basement/sco/morel/github/datasets/ensembl_8880_15/fasta_files/"
newick_dir = "/hits/basement/sco/morel/github/phd_experiments/datasets/treerecs/ensembl_8880_supports/support_values/"
alignment_file = "/hits/basement/sco/morel/github/phd_experiments/datasets/treerecs/ensembl_8880_supports/alignments.txt"
newick = "/hits/basement/sco/morel/github/phd_experiments/datasets/treerecs/ensembl_8880_supports/support_values/geneTrees.newick"

fastas = os.listdir(fasta_dir)

with open(alignment_file, "w") as f:
  f.write("GTR\n")
  for fasta in fastas:
    base = fasta.split(".")[0]



