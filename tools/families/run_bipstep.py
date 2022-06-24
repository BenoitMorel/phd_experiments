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
        writer.write(tree.write().replace("e-", "").replace("e+",""))
        writer.write("\n")
  return filepath


def exec_bipstep(gene_trees_file, output_dir, additional_arguments):
  command = []
  command.append(exp.bipstep_exec)
  command.append("-i")
  command.append(gene_trees_file)
  command.append("-p")
  prefix = os.path.join(output_dir, "bipstep")
  command.append(prefix)
  command.extend(additional_arguments) 
  print(" ".join(command))
  subprocess.check_call(command)
  return prefix + ".newick"

def run_bipstep(datadir, gene_trees, subst_model, additional_arguments):
  run_name = "bipstep"
  if ("--random-insertion" in additional_arguments):
    run_name += "rand"
  run_name += "_"
  run_name += gene_trees
  
  output_dir = fam.get_run_dir(datadir, subst_model, run_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  gene_trees_file = init_gene_trees_file(datadir, gene_trees, subst_model, output_dir)
  start = time.time()
  bipstep_tree = exec_bipstep(gene_trees_file, output_dir, additional_arguments)
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1, "runtimes") 
  shutil.copy(bipstep_tree, fam.get_species_tree(datadir, subst_model, run_name))

if (__name__ == "__main__"):
  if (len(sys.argv) < 4):
    print("Syntax python " + os.path.basename(__file__) + " datadir gene_trees subst_model")
    sys.exit(1)
  datadir = sys.argv[1]
  gene_trees = sys.argv[2]
  subst_model = sys.argv[3]
  additional_arguments = sys.argv[4:]
  run_bipstep(datadir, gene_trees, subst_model, additional_arguments)
  species_analyze.analyze(sys.argv[1])
  



