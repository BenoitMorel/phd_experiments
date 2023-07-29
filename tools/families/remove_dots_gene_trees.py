import os
import sys
import shutil
sys.path.insert(0, 'scripts')
import experiments as exp
import fam


def rename(datadir):
  for family in fam.get_families_list(datadir):
    p = fam.get_gene_tree_dir(datadir, family)
    for f in os.listdir(p):
      if (f.startswith("true")):
          continue
      sp = f.split(".")
      f2 = "_".join(sp[0:-2]) + "." + ".".join(sp[-2:])
      shutil.move(os.path.join(p,f), os.path.join(p,f2)) 
    #shutil.copyfile(s, d)
    #shutil.move(s, d)

if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " datadir")
    sys.exit(1)
  datadir = sys.argv[1]
  rename(datadir)



