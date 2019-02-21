import sys
import os
sys.path.insert(0, 'scripts')
import experiments as exp
import subprocess
import shutil

def run_stag(stag_map, stag_gene_tree_dir, output_species_tree):
  command = []
  command.append("python")
  command.append(exp.stag_script)
  command.append(stag_map)
  command.append(stag_gene_tree_dir)
  logs = subprocess.check_output(command)
  for line in logs.split("\n"):
    if (line.startswith("STAG species tree")):
      old_species_tree = line.split(" ")[3].replace("\n", "")
      shutil.copyfile(old_species_tree, output_species_tree)
      return
  print("Error: could not find stag output tree")

if __name__ == "__main__":
  if (len(sys.argv) != 4):
    print("syntax: python stag.py stag_map stag_gene_tree_dir output_species_tree")
    sys.exit(1)
  stag_map = sys.argv[1]
  stag_gene_tree_dir = sys.argv[2]
  output_species_tree = sys.argv[3]
  run_stag(stag_map, stag_gene_tree_dir, output_species_tree)



