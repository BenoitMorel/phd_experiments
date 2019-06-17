import os
import sys
sys.path.insert(0, os.path.join("tools", "trees"))
import create_random_tree
import fam


def run_random(datadir):
  for family in fam.get_families_list(datadir):
    msa = fam.get_alignment(datadir, family)
    random = fam.get_random_tree(datadir, family)
    create_random_tree.create_random_tree(msa, random)

if (__name__== "__main__"):
  max_args_number = 2
  if len(sys.argv) != max_args_number:
    print("Syntax error: python run_ALE.py datadir is_dna cores.")
    sys.exit(1)

  datadir = sys.argv[1]
  run_random(datadir)




