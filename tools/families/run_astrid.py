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
import species_analyze

def init_gene_trees_file(datadir, method, subst_model, output_dir):
  filepath = os.path.join(output_dir, "gene_trees.txt")
  with open(filepath, "w") as writer:
    for family in fam.get_families_list(datadir):
      raxml_tree_path = fam.build_gene_tree_path(datadir, subst_model, family, method)
      writer.write(open(raxml_tree_path).read().replace("\n", "") + "\n")
  return filepath

def init_mapping_file(datadir, output_dir):
  res = os.path.join(output_dir, "astrid_mappings.txt")
  with open(res, "w") as writer:
    dico = get_dico.get_species_to_genes(datadir)
    for species in dico:
      writer.write(species)
      writer.write(":")
      writer.write(",".join(dico[species]))
      writer.write("\n")
  return res

def exec_astrid(gene_trees_file, mapping_file, output_species_tree_file, mode):
  command = []
  command.append(exp.astrid_exec)
  command.append("-i")
  command.append(gene_trees_file)
  command.append("-a")
  command.append(mapping_file)
  command.append("-o")
  command.append(output_species_tree_file)
  if (mode == "bionj"):
    command.append("--bionj")
  elif (mode == "fastme"):
    command.append("-s")
  elif(mode == "default"):
    pass
  FNULL = open(os.devnull, 'w')
  res = subprocess.check_output(command)#, stderr=FNULL)

def run_astrid(datadir, method, subst_model, mode):
  name = "astrid-" + mode + "_" + method
  output_dir = fam.get_run_dir(datadir, subst_model, name + "_run")
  temp_species_tree = os.path.join(output_dir, "species_tree.newick")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  gene_trees_file = init_gene_trees_file(datadir, method, subst_model, output_dir)
  mapping_file = init_mapping_file(datadir, output_dir)
  print("Start executing astrid")
  start = time.time()
  exec_astrid(gene_trees_file, mapping_file, temp_species_tree, mode)
  output_species_tree = fam.get_species_tree(datadir, subst_model, name)
  shutil.move(temp_species_tree, output_species_tree)
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name(name, subst_model), time1, "runtimes") 

if (__name__ == "__main__"):
  if (len(sys.argv) < 4):
    print("Syntax python run_astrid.py datadir gene_trees subst_model [mode=default,bionj,fastme]")
    sys.exit(1)
  datadir = sys.argv[1]
  method = sys.argv[2]
  subst_model = sys.argv[3]
  mode = "default"
  if (len(sys.argv) == 5):
    mode = sys.argv[4]
  run_astrid(datadir, method, subst_model, mode)
  species_analyze.analyze(datadir) 
  




