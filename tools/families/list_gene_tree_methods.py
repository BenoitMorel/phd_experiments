import sys
import os
import fam

  
  
  
def list_gene_tree_methods(datadir):
  family = fam.get_families_list(datadir)[0]
  for f in os.listdir(fam.get_gene_tree_dir(datadir, family)):
    print(f)

if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " datadir")
    sys.exit(1)
  datadir = sys.argv[1]
  list_gene_tree_methods(datadir)



