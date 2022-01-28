import os
import sys
import shutil
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, 'scripts')
import experiments as exp
import fam
import contract_tree_branches
import numpy as np

def floatstr(my_float):
  return np.format_float_positional(my_float, trim='-')

def contract(datadir, gene_tree_method, subst_model, min_bl, min_support):
  method = gene_tree_method
  
  new_method = gene_tree_method 
  if (min_bl > 0.0):
    new_method = new_method + "-minbl" + floatstr(min_bl)
  if (min_support > 0.0):
    new_method = new_method + "-minsupport" + floatstr(min_support)
  for family in fam.get_families_list(datadir):
    src = fam.build_gene_tree_path(datadir, subst_model, family, method)
    dest = fam.build_gene_tree_path(datadir, subst_model, family, new_method)
    contract_tree_branches.contract_branches(src, dest, min_bl, min_support)    

if (__name__ == "__main__"):
  if (len(sys.argv) != 6):
    print("Syntax python " + os.path.basename(__file__) + " datadir gene_tree_method subst_model min_bl min_support")
    sys.exit(1)
  datadir = sys.argv[1]
  gene_tree_method = sys.argv[2]
  subst_model = sys.argv[3]
  min_bl = float(sys.argv[4])
  min_support = float(sys.argv[5])
  contract(datadir, gene_tree_method, subst_model, min_bl, min_support)



