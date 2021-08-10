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
import ete3

def init_gene_tree_dir(datadir, gene_trees, subst_model, output_dir):
  gene_tree_dir = os.path.join(output_dir, "gene_trees.txt")
  os.mkdir(gene_tree_dir)
  for family in fam.get_families_list(datadir):
    new_gene_tree_path = os.path.join(gene_tree_dir, family)
    with open(new_gene_tree_path, "w") as writer:
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
          species = m[leaf.name]
          leaf.name = leaf.name.replace("_", "DASHDASH")
          leaf.name = species + "_" + leaf.name
        writer.write(tree.write(format = 9))
  return gene_tree_dir


def exec_stride(species_tree_path, gene_trees_dir, output_dir, out_species_tree):
  command = []
  command.append(exp.python())
  command.append(exp.stride_script)
  command.append("-s")
  command.append("dash")
  command.append("-S")
  command.append(species_tree_path)
  command.append("-o")
  command.append(output_dir)
  command.append("-d")
  command.append(gene_trees_dir)
  FNULL = open(os.devnull, 'w')
  print(" ".join(command))
  subprocess.check_call(command, stderr=FNULL)
  species_tree = ete3.Tree(os.path.join(output_dir, "Species_tree_labelled.tre"), format=1)
  for leaf in species_tree:
    leaf.name = leaf.name.replace("DASHDASH", "_")
  species_tree.write(format=1, outfile = out_species_tree)

def run_stride(datadir, species_tree, gene_trees, subst_model):
  species_tree += "_" + gene_trees
  run_name = "stride_" + species_tree
  output_dir = fam.get_run_dir(datadir, subst_model, run_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  species_tree_path = fam.get_species_tree(datadir, subst_model, species_tree)
  gene_tree_dir = init_gene_tree_dir(datadir, gene_trees, subst_model, output_dir)
  start = time.time()
  out_species_tree = fam.get_species_tree(datadir, subst_model, run_name)
  exec_stride(species_tree_path, gene_tree_dir, output_dir, out_species_tree)
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1, "runtimes") 

if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax python run_stride.py datadir species_tree gene_trees subst_model")
    sys.exit(1)
  run_stride(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
  



