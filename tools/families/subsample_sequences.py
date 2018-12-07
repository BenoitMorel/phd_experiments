import os
import sys
import shutil
import ete3

def removeIfExists(directory):
  if (os.path.isdir(directory)):
    shutil.rmtree(directory)

def subsample_file(file_path, ratio):
  seqs_str = open(file_path).read()
  seqs = ete3.SeqGroup(seqs_str, format="fasta")
  for entry in seqs.get_entries():
    s = entry[1]
    s = s[:int(float(len(s)) * ratio)]
    seqs.set_seq(entry[0], s)
  seqs.write("fasta", file_path)

def subsample_sequences(input_dir, output_dir, ratio):
  if (os.path.isdir(output_dir)):
    print("Error, " + output_dir + " already exists")
    sys.exit(1)

  print("Copying files")
  shutil.copytree(input_dir, output_dir)
  removeIfExists(os.path.join(output_dir, "pargenes"))
  removeIfExists(os.path.join(output_dir, "treerecs_run"))
  removeIfExists(os.path.join(output_dir, "phyldog_run"))
  print("Subsampling")

  new_alignments_dir = os.path.join(output_dir, "alignments")
  for alignment in os.listdir(new_alignments_dir):
    subsample_file(os.path.join(new_alignments_dir, alignment), ratio)
  new_families_dir = os.path.join(output_dir, "families")
  for family in os.listdir(new_families_dir):
    subsample_file(os.path.join(new_families_dir, family, "alignment.msa"), ratio)

if (len(sys.argv) != 4):
  print("Syntax: python subsample_sequences.py dataset_dir new_dataset_dir ratio")
  sys.exit(1)

input_dir = sys.argv[1]
output_dir = sys.argv[2]
ratio = float(sys.argv[3])

subsample_sequences(input_dir, output_dir, ratio)
