import os
import sys
import shutil
sys.path.insert(0, 'scripts')
import experiments as exp
import fam

def remove_too_small_families(datadir):
  for family in fam.get_families_list(datadir):
    if (fam.get_family_genes_number(datadir, family) < 4):
      p = fam.get_family_path(datadir, family)
      shutil.rmtree(p)

if (__name__== "__main__"):
  max_args_number = 2
  if len(sys.argv) < max_args_number:
    print("Syntax error: python remove_too_small_families.py datadir.")
    sys.exit(1)
  remove_too_small_families(sys.argv[1])




