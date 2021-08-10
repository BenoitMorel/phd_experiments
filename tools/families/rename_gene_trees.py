import os
import sys
import shutil
sys.path.insert(0, 'scripts')
import experiments as exp
import fam


def rename(datadir, src, dest):
  for family in fam.get_families_list(datadir):
    p = fam.get_gene_tree_dir(datadir, family)
    s = os.path.join(p, src)
    d = os.path.join(p, dest)
    os.rename(s, d)

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python " + os.path.basename(__file__) + " datadir src_full_name dest_full_name")
    sys.exit(1)
  datadir = sys.argv[1]
  src = sys.argv[2]
  dest = sys.argv[3]
  rename(datadir, src, dest)


