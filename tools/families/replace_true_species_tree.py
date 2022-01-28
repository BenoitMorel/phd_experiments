import os
import sys
sys.path.insert(0, 'tools/families')
import fam
import shutil



if (__name__ == "__main__"): 
  if (len(sys.argv) < 2): 
    print("Syntax: python " + os.path.basename(__file__) + " species_tree datadirs")
    exit(1)
  species_tree = sys.argv[1]
  datadirs = sys.argv[2:]
  for datadir in datadirs:
    try:
      src = os.path.join(fam.get_species_dir(datadir), species_tree)
      dest = fam.get_true_species_tree(datadir)
      shutil.copyfile(src, dest)
    except:
      print("Failed to treat " + datadir)

