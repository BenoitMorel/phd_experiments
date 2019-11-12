import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/trees')
sys.path.insert(0, 'tools/mappings')
import experiments as exp
import fam
from read_tree import read_tree
import tree_format_converter
import get_dico
from ete3 import Tree

def get_new_species_name(old_species_name):
  return "ZXY" + old_species_name + "ZXY"

def get_old_species_name(new_species_name):
  return new_species_name.replace("ZXY", "")

def rename_gene_leaves(input_gene_tree, output_gene_tree, gene_to_species):
  tree = read_tree(input_gene_tree)
  species = {}
  for leaf in tree:
    species[gene_to_species[leaf.name]] = True
    leaf.name = get_new_species_name(gene_to_species[leaf.name]) + "_" + leaf.name
  if (len(species) < 4):
    return False
  open(output_gene_tree, "w").write(tree.write())
  return True

def init_gene_trees_file(datadir, subst_model, output_dir):
  gene_trees_file = os.path.join(output_dir, "gene_trees.txt")
  with open(gene_trees_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      newick_tree = fam.get_raxml_tree(datadir, subst_model, family)
      temp_newick_tree = "tempgenetree_" + family + ".newick"
      gene_to_species = get_dico.get_gene_to_species(datadir, family)
      ok = rename_gene_leaves(newick_tree, temp_newick_tree, gene_to_species)
      if (not ok):
        print("skipping family " + family)
        continue
      nexus_tree = "gene_trees_" + family + ".nexus"
      tree_format_converter.newick_to_nexus(temp_newick_tree, nexus_tree)
      writer.write(nexus_tree + "\n")
  return gene_trees_file

def init_species_file(datadir, output_dir):
  species_file = os.path.join(output_dir, "species.txt")
  leaves = read_tree(fam.get_species_tree(datadir)).get_leaves()
  leaves_str = []
  for leaf in leaves:
    leaves_str.append(get_new_species_name(leaf.name))
  with open(species_file, "w") as writer:
    writer.write("\n".join(leaves_str))
  return species_file

def build_guenomu_config_file(gene_trees_file, species_file, output_dir):
  config_file = os.path.join(output_dir, "config.txt")
  with open(config_file, "w") as writer:
    #writer.write("param_execute_action = import\n")
    writer.write("param_file_with_tree_files = " + gene_trees_file + "\n")
    writer.write("param_file_with_species_names = " + species_file + "\n")
    writer.write("param_reconciliation_prior = 0.0001\n")
    writer.write("param_use_distances = 1111000\n")
    writer.write("param_n_generations = 5000 10000\n")
    writer.write("param_n_samples = 100\n")
    writer.write("param_n_output = 2\n")
    writer.write("param_n_mc_samples = 10 \n")
    writer.write("param_mini_sampler = 0.8 5 5\n")
    writer.write("param_anneal = 10 1000 0.5 10\n")
    writer.write("\n")
    writer.write("\n")
  return config_file

def execute_guenomu(datadir, output_dir, config_file, step, cores):
  command = []
  
  command.append("mpiexec")
  command.append("-np")
  command.append(str(cores))
  command.append(exp.guenomu_exec)
  command.append("-z")
  command.append(str(step))
  command.append(config_file)
  subprocess.check_call(command)


def extract_results(datadir, subst_model):
  
  lines = open("job0.species.trprobs").readlines()
  translate = {}
  read_translate = False
  for line in lines:
    if ("Translate" in line):
      read_translate = True
      continue
    if (line.startswith(";")):
        read_translate = False
        continue
    if (read_translate):
      split = line.replace("\t", "").split(" ")
      translate[split[0]] = split[2].replace("\n", "").replace(",", "")
    if (line.startswith("tree tree_")):
      species_tree_str = line.split(" ")[-1].replace("\n", "")
      tree = Tree(species_tree_str)
      for leaf in tree:
        print(leaf.name)
        leaf.name = get_old_species_name(translate[leaf.name])
        print(leaf.name)
      open(fam.get_species_tree(datadir, subst_model, "guenomu"), "w").write(tree.write())
      print("saving in " + fam.get_species_tree(datadir, subst_model, "guenomu"))
      return



def run_guenomu(datadir, subst_model, cores):
  output_dir = fam.get_run_dir(datadir, subst_model, "guenomu_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  cwd = os.getcwd()
  try:
    os.chdir(output_dir)
    output_species_tree = fam.get_species_tree(datadir, subst_model, "guenomu")
    gene_trees_file = init_gene_trees_file(datadir, subst_model, output_dir)
    species_file = init_species_file(datadir, output_dir)
    config_file = build_guenomu_config_file(gene_trees_file, species_file, output_dir)
    execute_guenomu(datadir, output_dir, config_file, 0, cores)
    execute_guenomu(datadir, output_dir,  config_file, 1, cores)
    extract_results(datadir, subst_model)
    print("results in " + output_dir)
  finally:
    os.chdir(cwd)

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python run_guenomu.py datadir subst_model cores")
    sys.exit(1)
  run_guenomu(os.path.abspath(sys.argv[1]), sys.argv[2], int(sys.argv[3]))
  

