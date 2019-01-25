import sys
import os
import ete3

def phy_to_fasta_file(input_phy, output_fasta):
  lines = open(input_phy).readlines()
  phy_str = ""
  for line in lines:
    if ("_" in line):
      phy_str += (line[:-1] + " \n")
    else:
      phy_str += line 
  msa = ete3.SeqGroup(phy_str, format="phylip_relaxed")
  msa.write("fasta", output_fasta)

def phy_to_fasta_dir(input_phy_dir, output_fasta_dir):
  try:
    os.makedirs(output_fasta_dir)
  except:
    print("WARNING: directory " + output_fasta_dir + " already exists")
  for input_phy_base in os.listdir(input_phy_dir):
    input_phy = os.path.join(input_phy_dir, input_phy_base)
    output_fasta_base = ".".join(input_phy_base.split(".")[:-1]) + ".fasta"
    output_fasta = os.path.join(output_fasta_dir, output_fasta_base)
    phy_to_fasta_file(input_phy, output_fasta)


if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax: python phy_to_fasta.py input_phy_dir output_fasta_dir")
    exit(1)
  phy_to_fasta_dir(sys.argv[1], sys.argv[2])
