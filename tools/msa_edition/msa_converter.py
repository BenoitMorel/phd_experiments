import sys
import os
import ete3

def msa_convert(input_file, output_file, input_format, output_format, prefixes_dictionnary = None):
  lines = open(input_file).readlines()
  input_str = ""
  for line in lines:
    if ("_" in line):
      input_str += (line[:-1] + " \n")
    else:
      input_str += line 
  msa = ete3.SeqGroup(input_str, format=input_format) 
  if (prefixes_dictionnary != None):
    new_msa = ete3.SeqGroup()
    for seq in msa.get_entries():
      new_msa.set_seq(prefixes_dictionnary[seq[0]] + "_" + seq[0], seq[1])
    msa = new_msa 


  msa.write(output_format, output_file)

if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax: python phy_to_fasta.py file1 file2 format1 format2")
    print("Formats: phylip, phylip_relaxed, fasta, and all supported by ete3")
    exit(1)
  msa_convert(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])



