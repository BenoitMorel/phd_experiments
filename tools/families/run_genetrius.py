import os
import sys
import subprocess
import shutil
import time
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/trees')
sys.path.insert(0, 'tools/mappings')
import experiments as exp
import fam
from read_tree import read_tree
import saved_metrics
import get_dico
import ete3 
  


def get_matrix(datadir, species_tree):
  matrix_path = os.path.join(datadir, "misc", "presence_matrix.txt")
  labels = ete3.Tree(species_tree).get_leaf_names()
  N = len(labels)
  labelToId = {}
  for spid in range(0, N):
    labelToId[labels[spid]] = spid
  families = fam.get_families_list(datadir)
  K = len(families)
  matrix = []
  for spid in range(0, N):
    matrix.append(["0"] * K)
  for k in range(0, K):
    family = families[k]
    d = get_dico.get_species_to_genes_family(datadir, family)
    for species in d:
      spid = labelToId[species]
      matrix[spid][k] = "1"

  with open(matrix_path, "w") as writer:
    writer.write(str(N))
    writer.write(" ")
    writer.write(str(K))
    writer.write("\n")
    for spid in range(0, N):
      writer.write(str(labels[spid]))
      writer.write(" ")
      writer.write(" ".join(matrix[spid]))
      writer.write("\n")
  return matrix_path
  
  
def execute(species_tree, matrix, cores):
  command = []
  command.append(exp.genetrius_exec)
  command.append("-nt")
  command.append(str(cores))
  command.append("--gentrius")
  command.append(species_tree)
  command.append("-pr_ab_matrix")
  command.append(matrix)
  print(" ".join(command))
  subprocess.check_call(command)

def run_genetrius(datadir, species_tree, cores):
  matrix = get_matrix(datadir, species_tree)
  execute(species_tree, matrix, cores)

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python " + os.path.basename(__file__) + " species_tree cores")
    sys.exit(1)
  species_tree= sys.argv[1]
  cores = sys.argv[2]
  datadir = species_tree.split("/species_trees/")[0]
  run_genetrius(datadir, species_tree, cores)

