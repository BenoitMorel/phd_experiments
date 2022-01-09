import os
import sys
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "mappings"))
import create_random_tree
import fam
import get_dico

def run_random(datadir):
  for family in fam.get_families_list(datadir):
    genes = get_dico.get_genes(datadir, family) 
    random = fam.get_random_tree(datadir, family)
    tree = create_random_tree.create_random_tree_from_species(genes)
    tree.write(outfile = random)

if (__name__== "__main__"):
  max_args_number = 2
  if len(sys.argv) != max_args_number:
    print("Syntax error: python run_random.py datadir.")
    sys.exit(1)

  datadir = sys.argv[1]
  run_random(datadir)




