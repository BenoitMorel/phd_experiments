import os
import sys
sys.path.insert(0, os.path.join("scripts"))
sys.path.insert(0, os.path.join("tools", "families"))
import fam
import experiments as exp
import shutil
import extract



def extract_mammals():
  src = fam.get_datadir("tuto_mammals")
  dest = os.path.join(exp.workshop_root, "data", "practical3", "mammals")
  families = fam.get_families_list(src)
  print(len(families))
  # I generated trento_mammals with the script current_experiments/trento/generate_mammals.sh
  # this step only extracts it to the workshop format
  extract.extract(src, dest, families, extract_mappings = True)

if (__name__ == "__main__"): 
  extract_mammals()


