import sys
import os
import ete3

def write_nexus(msa, output):
  with open(output, "w") as writer:
    ntax = len(msa.get_entries())
    first_seq = msa.get_entries()[0][1]
    nchar = len(first_seq)
    datatype = "dna"
    dna_chars = set(["a", "c", "g", "t"])
    for c in first_seq.lower():
      if (not c in dna_chars):
        datatype = "protein"
        break
    writer.write("#NEXUS\n")
    writer.write("begin data;\n")
    writer.write("\tdimensions ntax=" + str(ntax) + " nchar= " + str(nchar) + ";\n")
    writer.write("\tformat datatype=" + datatype + " interleave=no gap=-;\n")
    writer.write("\tmatrix\n")
    for seq in msa.get_entries():
      writer.write(seq[0] + "\t" + seq[1] + "\n") 
    writer.write("\t;\n")
    writer.write("end;")

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

  if (output_format.lower() == "nexus"):
    write_nexus(msa, output_file)
  else:
    msa.write(output_format, output_file)

if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax: python phy_to_fasta.py file1 file2 format1 format2")
    print("Formats: phylip, phylip_relaxed, fasta, and all supported by ete3")
    exit(1)
  msa_convert(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])



