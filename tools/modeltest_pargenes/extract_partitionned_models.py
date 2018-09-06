import argparse
import sys
import os

# useful to get partition files with modeltest models from
# a pargenes run on a set of phy files obtained with
# the script tools/msa_edition/split_phy_per_partitions.py

parser = argparse.ArgumentParser()
parser.add_argument('-p', "--partition-file", 
    dest="partition_file", 
    help="Partitions file")
parser.add_argument('-d', "--pargenes-dir", 
    dest="pargenes_dir", 
    help="ParGenes run directory")
parser.add_argument('-o', "--output", 
    dest="output", 
    help="Output partition file")
op = parser.parse_args()

def get_model_from_log(log_file, modeltest_criteria):
  with open(log_file) as reader:
    read_next_model = False
    for line in reader.readlines():
      if (line.startswith("Best model according to " + modeltest_criteria)):
          read_next_model = True
      if (read_next_model and line.startswith("Model")):
        model = line.split(" ")[-1][:-1]
        return model
  return None


partition_lines = open(op.partition_file).readlines()
modeltest_results = os.path.join(op.pargenes_dir, "modeltest_run", "results")
criteria = ["AIC", "AICc", "BIC"]

for criterion in criteria:
  writer = open(op.output + "." + criterion, "w")
  for line in partition_lines:
    split = line.split(" ")
    name = split[1] + "_phy"
    modeltest_outfile = os.path.join(modeltest_results, name, name + ".out")  
    model = get_model_from_log(modeltest_outfile, criterion)
    split[0] = model + ","
    writer.write(" ".join(split))



