import sys
import os
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import experiments as exp
import fam


class SeqEntry():
  def __init__(self, line, family):
    split = line.split(" ")
    temp = split[1].split("_")
    self.species = temp[0].title() + "_" + temp[1].title()
    self.gene = split[2]
    self.chromozome = split[3]
    self.begin = int(split[4])
    self.end = int(split[5])
    self.wtf1 = split[6]
    self.wtf2 = split[7]
    self.wtf_function = split[8][:-1]
    self.family = family

  def __str__(self):
    res = "("
    res += "species:" + self.species + ", "
    res += "gene:" + self.gene + ", "
    res += "chrom:" + self.chromozome + ", "
    res += "begin:" + str(self.begin) + ", "
    res += "end:" + str(self.end) + ", "
    res += "family:" + self.family
    res += ")"
    return res

def parse_nhx_emf(emf_file):
  per_family_seq_entries = {}  
  family_index = 0
  family = "family_" + str(family_index)
  per_family_seq_entries[family] = []
  for line in open(emf_file).readlines():
    if (line[:3] == "SEQ"):
      entry = SeqEntry(line, family)
      per_family_seq_entries[entry.family].append(entry)
    elif (line[:4] == "DATA"):
      family_index += 1
      family = "family_" + str(family_index)
      per_family_seq_entries[family] = []
  print("number of gene families: " + str(len(per_family_seq_entries)))
  return per_family_seq_entries


def parse_fasta(fasta_file):
  alignments_dico = {}
  cur_name = ""
  cur_seq = ""
  for line in open(fasta_file).readlines():
    if (line.strip() == "" or line[0] == "/"):
      continue
    if (line[0] == ">"):
      if (len(cur_name)):
        alignments_dico[cur_name] = cur_seq
      cur_name = line[1:-1]
      cur_seq = line
    else:
      cur_seq += line
  alignments_dico[cur_name] = cur_seq
  print("number of DNA seq: " + str(len(alignments_dico)))
  return alignments_dico

def export_msa(seq_entries, alignments_dico, output_file):
  with open(output_file, "w") as writer:
    for seq_entry in seq_entries:
      writer.write(alignments_dico[seq_entry.gene])


def export(per_family_seq_entries, alignments_dico, datadir):
  os.makedirs(datadir)
  families_dir = fam.getFamilies(datadir)
  os.makedirs(families_dir)
  for family in per_family_seq_entries:
    family_dir = fam.getFamily(datadir, family)
    os.makedirs(family_dir)
    seq_entries = per_family_seq_entries[family]
    export_msa(seq_entries, alignments_dico, fam.getAlignment(datadir, family)) 
    
def extract_from_ensembl(nhx_emf_file, fasta_file, datadir):
  per_family_seq_entries = parse_nhx_emf(nhx_emf_file)
  alignments_dico = parse_fasta(fasta_file)
  export(per_family_seq_entries, alignments_dico, datadir)

if (__name__ == "__main__"): 
  if (len(sys.argv) != 4): 
    print("Syntax: python " + os.path.basename(__file__) + " nhx_emf fasta output")
    exit(1)
  nhx_emf_file = sys.argv[1]
  fasta_file = sys.argv[2]
  datadir = sys.argv[3]
  extract_from_ensembl(nhx_emf_file, fasta_file, datadir)

