
import os
import sys
from ete3 import SeqGroup

class Partition:
  def __init__(self):
    self.start = 0
    self.width = 0
    self.name = ""

  def load_from_line(self, line):
    print(line)
    line = line.replace(",", "")
    line = line.replace("DNA", "")
    line = line.replace("charset", "")
    line = line.replace("=", " ")
    line = line.replace(" - ", "-")
    line = line.replace(";", "")
    split = line.split()
    print(split)
    self.name = split[0]
    bounds = split[1].split("-")
    self.start = int(bounds[0]) - 1
    self.width = int(bounds[1]) - self.start
    print(str(self.start) + " " + str(self.width))

def load_partitions(file_name):
  partitions = []
  lines = open(file_name).readlines()
  for line in lines:
    if (not "-" in line):
      continue
    part = Partition()
    part.load_from_line(line)
    partitions.append(part)
  return partitions


def read_alignment(alignment_file):
  msa = None
  formats = ["fasta", "phylip_relaxed", "iphylip_relaxed", "phylip_interleaved"]
  for f in formats:
    try:
      msa = SeqGroup(alignment_file, f)
      return msa
    except:
      pass

def get_sub_alignment(partition, concatenated_alignment):
  sub_alignment = SeqGroup()
  
  for entry in concatenated_alignment.get_entries():
    seq = entry[1]
    sub_seq = seq[partition.start:partition.start + partition.width]
    sub_alignment.set_seq(entry[0], sub_seq)
  return sub_alignment
 
def get_sub_alignments(alignment_file, partition_file):
  concatenated_alignment = read_alignment(alignment_file)
  partitions = load_partitions(partition_file)
  res = {}
  for partition in partitions:
    sub_alignment = get_sub_alignment(partition, concatenated_alignment)
    res[partition.name] = sub_alignment
  return res

def split_alignment(alignment_file, partition_file, output_dir):
  os.mkdir(output_dir)
  concatenated_alignment = read_alignment(alignment_file)
  partitions = load_partitions(partition_file)
  for partition in partitions:
    output_file = os.path.join(output_dir, partition.name + ".fasta")
    sub_alignment = get_sub_alignment(partition, concatenated_alignment)
    sub_alignment.write("fasta", output_file)

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python " + os.path.basename(__file__) + " msa partition output_msa")
    sys.exit(1)
  split_alignment(sys.argv[1], sys.argv[2], sys.argv[3])
