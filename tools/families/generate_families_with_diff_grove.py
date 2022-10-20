import sys
import os
import shutil
import random
import glob
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/database')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
from find_diff_datasets_julia import get_diff_and_propmax
from ete3 import SeqGroup



def get_nt_from_logs(logfile):
  lines = open(logfile).readlines()
  for line in lines:
    if ("Loaded alignment with" in line):
      sp = line.split()
      print(sp[4] + " " +sp[7])
      return int(sp[4]) * int(sp[7])
  return 0


def generate(inputdir, min_var, max_propmax, max_fam, max_nt, seed):
  random.seed(seed)
 
  outputdir = "diffgrove" 
  outputdir += "_var" + str(min_var)
  outputdir += "_propmax" + str(max_propmax)
  outputdir += "_fam" + str(max_fam)
  outputdir += "_nt"
  outputdir += "_seed" + str(seed)
  outputdir = fam.get_datadir(outputdir)
  fam.init_top_directories(outputdir)

  all_families = os.listdir(inputdir)
  all_families.sort()
  random.shuffle(all_families)
  print("Total number of families " + str(len(all_families))) 
  families = []
  index = 0
  for family in all_families:
    index += 1
    famdir = os.path.join(inputdir, family)
    #if (not os.path.isfile(os.path.join(famdir, "data.sqlite3"))):
    #  continue
    raxmldir = os.path.join(inputdir, family, "output_files", "raxmlng", "evaluation")
    
    pat = os.path.join(raxmldir, "*raxml.log")
    print(pat)
    anylog = None
    try:
      anylog = glob.glob(pat)[0]
    except:
      continue
    bases = get_nt_from_logs(anylog)
    print(bases)
    if (bases > max_nt):
      continue

    t = None
    try:
      t = get_diff_and_propmax(famdir)
    except:
      continue
    print(t)
    var = t[0]
    propmax = t[1]
    
    if (var >= min_var and propmax <= max_propmax):
      families.append(family)
      print("*********************")
      print("Add " + family + " " + str(len(families)))
      print(" after " + str(index) + " iterations")
    if (len(families) >= max_fam):
      break

  for family in families:
    raxmldir = os.path.join(inputdir, family, "output_files", "raxmlng", "evaluation")
   
    ali = None
    try:
      pat = os.path.join(raxmldir, "*.phy")
      ali =   glob.glob(pat)[0]
    except:
      continue
    fam.init_family_directories(outputdir, family)
    shutil.copy(ali, fam.get_alignment(outputdir, family))
    trees = glob.glob(os.path.join(raxmldir, "*.bestTree"))
    print(len(trees))
    with open(fam.get_raxml_multiple_trees(outputdir, "GTR+G", family, 100), "w") as writer:
      for tree in trees:
        writer.write(open(tree).read())

  #fam.postprocess_datadir(outputdir)


if (__name__ == "__main__"): 
  if (len(sys.argv) < 6): 
    print("Syntax: python " + os.path.basename(__file__) + " inputdir min_var max_propmax max_families max_bases seed")
    exit(1)
  inputdir = sys.argv[1]
  min_var = float(sys.argv[2])
  max_propmax = float(sys.argv[3])
  max_fam = int(sys.argv[4])
  max_nt = int(sys.argv[5])
  seed = int(sys.argv[6])
  generate(inputdir, min_var, max_propmax, max_fam, max_nt, seed)

