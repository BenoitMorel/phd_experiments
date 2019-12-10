import sys
import os
import shutil
import functools
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
import create_random_tree



def generate_from_msas(msas_dir, datadir):
  fam.init_top_directories(datadir)
  
  families = []
  for f in os.listdir(msas_dir):
    families.append(f.split(".")[0])
  fam.init_families_directories(datadir, families)
  for f in os.listdir(msas_dir):
    src = os.path.join(msas_dir, f)
    dest = fam.get_alignment(datadir, f.split(".")[0])
    shutil.copyfile(src, dest)
  fam.postprocess_datadir(datadir)


if (__name__ == "__main__"): 
  if (len(sys.argv) < 3): 
    print("Syntax: python " + os.path.basename(__file__) + " msas_dir datadir")
    exit(1)
  msas_dir = sys.argv[1]
  datadir = sys.argv[2]
  generate_from_msas(msas_dir, datadir)

