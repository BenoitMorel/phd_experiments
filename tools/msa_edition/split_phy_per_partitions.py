import sys
import os

input_phylip = sys.argv[1]
input_part = sys.argv[2]
output_dir = sys.argv[3]

class Partition:
  def __init__(self):
    self.start = 0
    self.width = 0
    self.name = ""

  def load_from_line(self, line):
    split = line.split(" ")
    self.name = split[1]
    bounds = split[3].split("-")
    self.start = int(bounds[0]) - 1
    self.width = int(bounds[1][:-1]) - self.start
  

class Sequence:
  def __init__(self, label, sequence):
    self.label = label
    self.sequence = sequence

class MSA:
  def __init__(self, name):
    self.sequences = []
    self.width = 0
    self.name = name

  def add_sequence(self, sequence):
    self.sequences.append(sequence)

  def load_phylip_sequential(self, file_name):
    lines = open(file_name).readlines()
    split = lines[0].split(" ")
    taxa = int(split[0])
    self.width = int(split[1][:-1])
    for line in lines[1:]:
      split = line.split(" ")
      self.add_sequence(Sequence(split[0], split[1][:-1]))

  def save_phylip_sequential(self, file_name):
    with open(file_name, "w") as f:
      f.write(str(len(self.sequences)) + " " + str(self.width) + "\n")
      for seq in self.sequences:
        f.write(seq.label)
        f.write(" ")
        f.write(seq.sequence)
        f.write("\n")


def load_partitions(file_name):
  partitions = []
  lines = open(file_name).readlines()
  for line in lines:
    part = Partition()
    part.load_from_line(line)
    partitions.append(part)
  return partitions

def get_sub_seq(seq, partition):
  return Sequence(seq.label, seq.sequence[partition.start: partition.start + partition.width])

def get_sub_msa(msa, partition):
  submsa = MSA(partition.name)
  submsa.width = partition.width
  for seq in msa.sequences:
    submsa.add_sequence(get_sub_seq(seq, partition))
  return submsa

def split_per_partition(msa, partitions):
  msas = []
  for partition in partitions:
    msas.append(get_sub_msa(msa, partition))
  return msas

def save_msas(msas, output_dir):
  for msa in msas:
    file_name = os.path.join(output_dir, msa.name + ".phy")
    msa.save_phylip_sequential(file_name)
      

msa = MSA("MSA")
msa.load_phylip_sequential(input_phylip)

partitions = load_partitions(input_part)

submsas = split_per_partition(msa, partitions)

try:
  os.makedirs(output_dir)
except:
  pass

save_msas(submsas, output_dir)

