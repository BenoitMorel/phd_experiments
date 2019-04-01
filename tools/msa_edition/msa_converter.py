import sys
import os
import ete3
  
def msa_convert(input_file, output_file, input_format, output_format):
  lines = open(input_file).readlines()
  input_str = ""
  for line in lines:
    if ("_" in line):
      input_str += (line[:-1] + " \n")
    else:
      input_str += line 
  msa = ete3.SeqGroup(input_str, format=input_format)
  msa.write(output_format, output_file)

if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax: python phy_to_fasta.py file1 file2 format1 format2")
    print("Formats: phylip, phylip_relaxed, fasta, and all supported by ete3")
    exit(1)
  msa_convert(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

