import sys
import os
import shutil
import subprocess
import functools
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
import create_random_tree
from ete3 import Tree
from ete3 import SeqGroup
import re


def generate(gene_trees_dir, species_tree_path, dict_path, datadir):
  d = {}
  for line in open(dict_path).readlines():
    sp = line.replace("\n", "").split() 
    print(line)



if (__name__ == "__main__"): 
  if (len(sys.argv) != 5): 
    print("Syntax: python " + os.path.basename(__file__) + " gene_trees_dict species_tree dictionary output_dir")
    exit(1)
  gene_trees_dir = sys.argv[1]
  species_tree_path = sys.argv[2]
  dict_path = sys.argv[3]
  datadir = sys.argv[4]
  generate(gene_trees_dir, species_tree_path, dict_path, datadir)



