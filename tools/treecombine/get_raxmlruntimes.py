import os
import sys
sys.path.insert(0, 'tools/families')
import fam
import glob

def get_runtime_from_logs(logfile):
  lines = open(logfile).readlines()
  for line in lines:
    if ("Elapsed time" in line):
      return float(line.split(" ")[2])

def get_raxmlruntimes(datadir, databasedir):
  runtimes = []
  for family in fam.get_families_list(datadir):
    print(family)
    raxmldir = os.path.join(databasedir, family, "output_files", "raxmlng", "inference")
    pat = os.path.join(raxmldir, "*raxml.log")
    fam_runtime = 0
    for log in glob.glob(pat):
      fam_runtime += get_runtime_from_logs(log)
    print(fam_runtime)
    runtimes.append(fam_runtime)
  total_runtime = sum(runtimes)
  print("Total runtime = " + str(total_runtime) + "s")


if (__name__ == "__main__"):
  if (len(sys.argv) < 3):
    print("Syntax python " + os.path.basename(__file__) + " datadir treebase_dir ")
    sys.exit(1)
  datadir = sys.argv[1]
  databasedir = sys.argv[2]
  get_raxmlruntimes(datadir, databasedir)




