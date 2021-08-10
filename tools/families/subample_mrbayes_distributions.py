import os
import sys
import subprocess
import shutil
import time
import concurrent.futures
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "msa_edition"))
import saved_metrics
from run_mrbayes import MrbayesInstance
import experiments as exp

def substample_distribution(src, dest, reduce_by):
  input_lines = open(src).readlines()
  idx = 0
  with open(dest, "w") as writer:
    for line in input_lines:
      if (idx % reduce_by == 0):
        writer.write(line)
      idx += 1


def subsample(datadir, gene_trees, reduce_by):
  inst = MrbayesInstance.get_instance(datadir, gene_trees) 
  old_tag = inst.get_tag()
  inst.frequency = inst.frequency / reduce_by
  inst.burnin = inst.burnin / reduce_by
  subst_model = inst.subst_model
  new_tag = inst.get_tag()
  print(new_tag + " -> " + old_tag)
  for family in fam.get_families_list(datadir):
    
    gene_trees = fam.build_gene_tree_path(datadir, subst_model, family, old_tag)
    new_gene_trees = fam.build_gene_tree_path(datadir, subst_model, family, new_tag)
    substample_distribution(gene_trees, new_gene_trees, reduce_by) 


if (__name__== "__main__"):
  if len(sys.argv) < 4:
    print("Syntax error: python " + os.path.basename(__file__) + " gene_tree reduce_by (for instance 10) datadir_list")
    print(len(sys.argv))
    sys.exit(0)

  gene_trees = sys.argv[1]
  reduce_by = int(sys.argv[2])
  datadirs = sys.argv[3:]
  for datadir in datadirs:
    print(datadir)
    subsample(datadir, gene_trees, reduce_by)

