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
from read_tree import read_trees_list
import saved_metrics
import get_dico
import species_analyze

def init_gene_trees_file(datadir, gene_trees, subst_model, output_dir):
  filepath = os.path.join(output_dir, "gene_trees.txt")
  print(filepath) 
  with open(filepath, "w") as writer:
    for family in fam.get_families_list(datadir):
      gene_tree_path = fam.build_gene_tree_path(datadir, subst_model, family, gene_trees)
      mapping_path = fam.get_mappings(datadir, family)
      m = get_dico.get_gene_to_species(datadir, family)
      for line in open(mapping_path).readlines():
        split = line.replace("\n", "").split(":")
        species = split[0]
        genes = split[1].split(";")
        for gene in genes:
          m[gene] = species
      trees = read_trees_list(gene_tree_path)
      for tree in trees:
        for leaf in tree:
          leaf.name = m[leaf.name]
        writer.write(tree.write().replace("e-", ""))
        writer.write("\n")
  return filepath


def exec_astral(gene_trees_file, output_dir, addition_species_tree):
  command = []
  command.append("java")
  command.append("-Xms700G")
  command.append("-Xmx700G")
  command.append("-jar")
  command.append(exp.astralmp_jar)
  command.append("-i")
  command.append(gene_trees_file)
  command.append("-o")
  out = os.path.join(output_dir, "out.txt")
  command.append(os.path.join(output_dir, "out.txt"))
  if (addition_species_tree != None):
    command.append("-e")
    command.append(addition_species_tree)
  #command.append("-r")
  #command.append("all")
  print(" ".join(command))
  FNULL = open(os.devnull, 'w')
  res = subprocess.check_output(command)
  #res = subprocess.check_output(command, stderr=FNULL)
  print(res)
  return out

def run_astral(datadir, gene_trees, subst_model, addition_species_tree = None):
  run_name = "astralmp_" + gene_trees
  if (addition_species_tree != None):
    run_name += "-additional"
  output_dir = fam.get_run_dir(datadir, subst_model, run_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  gene_trees_file = init_gene_trees_file(datadir, gene_trees, subst_model, output_dir)
  start = time.time()
  astral_tree = exec_astral(gene_trees_file, output_dir, addition_species_tree)
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1, "runtimes") 
  shutil.copy(astral_tree, fam.get_species_tree(datadir, subst_model, run_name))

if (__name__ == "__main__"):
  if (len(sys.argv) < 4):
    print("Syntax python " + os.path.basename(__file__) + " datadir gene_trees subst_model")
    sys.exit(1)
  datadir = sys.argv[1]
  gene_trees = sys.argv[2]
  subst_model = sys.argv[3]
  addition_species_tree = None
  if (len(sys.argv) > 4):
      addition_species_tree = sys.argv[4]
  run_astral(datadir, gene_trees, subst_model, addition_species_tree)
  species_analyze.analyze(sys.argv[1])
  


