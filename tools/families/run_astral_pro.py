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
def init_gene_trees_file(datadir, method, subst_model, output_dir):
  filepath = os.path.join(output_dir, "gene_trees.txt")
  with open(filepath, "w") as writer:
    for family in fam.get_families_list(datadir):
      gene_tree_path = fam.build_gene_tree_path(datadir, subst_model, family, method)
      tree = ete3.Tree(gene_tree_path)
      for node in tree.traverse("postorder"):
        node.support = 100
        node.dist = 1.0
      towrite = tree.write()
      while (towrite[-1] == "\n"):
        towrite = towrite[:-1] 
      writer.write(towrite + "\n")
  return filepath

def init_mapping_file(datadir, output_dir):
  res = os.path.join(output_dir, "astralpro_mappings.txt")
  with open(res, "w") as writer:
    dico = get_dico.get_species_to_genes(datadir)
    for species in dico:
      for gene in dico[species]:
        writer.write(gene + " " + species + "\n")
  return res

def exec_astralpro(gene_trees_file, mapping_file, output_species_tree_file):
  command = []
  tmp_output_species_tree_file = output_species_tree_file + ".tmp"
  command.append("java")
  command.append("-Xms500G")
  command.append("-Xmx500G")
  command.append("-Djava.library.path=" + exp.astralpro_root + "/ASTRAL-MP/lib")
  command.append("-jar")
  command.append(exp.astralpro_jar)
  command.append("-i")
  command.append(gene_trees_file)
  command.append("-a")
  command.append(mapping_file)
  command.append("-o")
  command.append(tmp_output_species_tree_file)
  command.append("--seed")
  command.append("692")
  command.append("-t")
  command.append("1")
  FNULL = open(os.devnull, 'w')
  subprocess.check_call(command)
  shutil.move(tmp_output_species_tree_file, output_species_tree_file)
  

def run_astralpro(datadir, method, subst_model):
  print("Start astral pro script")
  run_name = "astralpro_" + method
  output_dir = fam.get_run_dir(datadir, subst_model, run_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  gene_trees_file = init_gene_trees_file(datadir, method, subst_model, output_dir)
  mapping_file = init_mapping_file(datadir, output_dir)
  print("Start executing astralpro")
  start = time.time()
  exec_astralpro(gene_trees_file, mapping_file, fam.get_species_tree(datadir, subst_model, run_name))
  time1 = (time.time() - start)
  print("Runtime: " + str(time1) + "s")
  saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1, "runtimes") 

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python " + os.path.basename(__file__) + " datadir gene_trees subst_model")
    sys.exit(1)
  run_astralpro(sys.argv[1], sys.argv[2], sys.argv[3])
  



