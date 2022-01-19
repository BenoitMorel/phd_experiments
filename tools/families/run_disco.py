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
import get_dico
from read_tree import read_trees_list
from ete3 import Tree

SEP = "|"

def exec_disco(src, dest):
  command = []
  command.append(exp.python3())
  command.append(exp.disco_script)
  command.append("-i")
  command.append(src)
  command.append("-o")
  command.append(dest)
  command.append("-d")
  command.append(SEP)
  res = subprocess.check_output(command)
  
def prefix_with_species(src, dest, dico):
  tree = Tree(src, 1)
  for leaf in tree.get_leaves():
    species = dico[leaf.name]
    leaf.name = species + SEP + leaf.name 
  tree.write(outfile = dest)

def remove_prefix(src, dest):
  total = 0
  with open(dest, "w") as writer:
    
    for line in open(src).readlines():
      if (len(line) < 3):
        continue
      total = total + 1
      tree = Tree(line, 1)
      for leaf in tree.get_leaves():
        sp = leaf.name.split(SEP)
        leaf.name = SEP.join(sp[1:])
      writer.write(tree.write())
      writer.write("\n")
  return total

def run_disco(datadir, gene_trees, subst_model):
  total = 0
  for family in fam.get_families_list(datadir):
    f1 = fam.build_gene_tree_path(datadir, subst_model, family, gene_trees)
    f2 = f1 + ".temp"
    f3 = f2 + ".temp"
    f4 = fam.build_gene_tree_path(datadir, subst_model, family, gene_trees + "-disco")
    dico = get_dico.get_gene_to_species(datadir, family)
    prefix_with_species(f1, f2, dico)
    exec_disco(f2, f3)
    total = total + remove_prefix(f3, f4)
    os.remove(f2)
    os.remove(f3)
    print("Number of trees: " + str(total))

if (__name__ == "__main__"):
  if (len(sys.argv) < 4):
    print("Syntax python " + os.path.basename(__file__) + " datadir gene_trees subst_model")
    sys.exit(1)
  datadir = sys.argv[1]
  gene_trees = sys.argv[2]
  subst_model = sys.argv[3]
  run_disco(datadir, gene_trees, subst_model)
  



