import os
import sys
import extract_gene_family as extract

def get_msa_id(msa_file):
  return msa_file.replace(".", "_")

if (len(sys.argv) != 3):
  print("Syntax : input_dataset_dir output_dir")
  sys.exit(1)

input_dir = sys.argv[1]
output_dir = sys.argv[2]

try:
  os.makedirs(output_dir)
except:
  print("Could not create directory " + output_dir)
  print("Please check that it does not exist yet")
  sys.exit(1)

input_msas_dir = os.path.join(input_dir, "alignments")

for msa in os.listdir(input_msas_dir):
  family_output_dir = os.path.join(output_dir, get_msa_id(msa))
  extract.extract_gene_family(input_dir, msa, family_output_dir)
  



