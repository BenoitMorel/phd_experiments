import shutil
import sys
import os
import subprocess
import time
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/mappings')
import get_dico
import experiments as exp
import fam
import species_analyze
import saved_metrics

def init_gene_trees_file(datadir, gene_trees, subst_model, output_dir):
  filepath = os.path.join(output_dir, "gene_trees.txt")
  with open(filepath, "w") as writer:
    for family in fam.get_families_list(datadir):
      gene_tree_path = fam.build_gene_tree_path(datadir, subst_model, family, gene_trees)
      towrite = open(gene_tree_path).read()
      while (towrite[-1] == "\n"):
        towrite = towrite[:-1]
      writer.write(towrite + "\n")
  return filepath

def init_mapping_file(datadir, output_dir):
  res = os.path.join(output_dir, "njrax_mappings.txt")
  with open(res, "w") as writer:
    dico = get_dico.get_species_to_genes(datadir)
    for species in dico:
      for gene in dico[species]:
        writer.write(gene + " " + species + "\n")
  return res

def exec_njrax(algo, gene_trees_file, mapping_file, output_tree):
  command = []
  command.append(exp.njrax_exec)
  command.append(algo)
  command.append(gene_trees_file)
  command.append(mapping_file)
  command.append(output_tree)
  subprocess.check_call(command)
    

def run_njrax(datadir, algo, gene_trees, subst_model):
  run_name = "njrax-" + algo + "-" + gene_trees
  output_dir = fam.get_run_dir(datadir, subst_model, run_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  gene_trees_file = init_gene_trees_file(datadir, gene_trees, subst_model, output_dir)
  output_tree = fam.get_species_tree(datadir, subst_model, run_name) 
  mapping_file = init_mapping_file(datadir, output_dir)
  start = time.time()
  exec_njrax(algo, gene_trees_file, mapping_file, output_tree)
  time1 = (time.time() - start)
  print("Runtime: " + str(time1) + "s")
  saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1, "runtimes") 

  


if (__name__== "__main__"):
  max_args_number = 5
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_njrax.py datadir algo gene_trees subst_models.")
    sys.exit(0)


  datadir = sys.argv[1]
  algo = sys.argv[2]
  gene_trees = sys.argv[3]
  subst_model = sys.argv[4] 
  run_njrax(datadir, algo, gene_trees, subst_model)
  species_analyze.analyze(datadir)


