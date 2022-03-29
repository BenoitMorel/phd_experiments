import sys
import os
from ete3 import Tree
sys.path.insert(0, os.path.join("tools", "trees"))
import grf_distance
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "print"))
sys.path.insert(0, os.path.join("script"))
import saved_metrics
import fam
from aligned_printer import AlignedPrinter
import experiments as exp

def get_runs(datadir):
  res = []
  for f in os.listdir(fam.get_species_dir(datadir)):
    split = f.split(".")
    run = ".".join(split[0:-3])
    if (len(split) >= 4):
      res.append(run + "." + split[-3])
  return res

def get_species_tree(datadir, run):
  return os.path.join(fam.get_species_dir(datadir), run + ".speciesTree.newick")


def analyze(datadir):
  runs = get_runs(datadir)
  true_tree = fam.get_species_tree(datadir)
  trees = {}
  for run in runs:
    trees[run] = get_species_tree(datadir, run)

  rf_printer = AlignedPrinter()
  grf_printer = AlignedPrinter()

  for run in trees:
    #try:
      tree = trees[run]
      if (tree == None):
        continue
      # compute RF distances
      (rf, grf) = grf_distance.compute_distances(true_tree, tree)
      # relativee distances 
      
      # save metrics
      saved_metrics.save_metrics(datadir, run, str(grf), "species_grf") 
      # add to printer
      rf_printer.add(run + ":", str(rf))
      grf_printer.add(run + ":", str(grf))

  print("")
  print("Unrooted average RF:")
  rf_printer.sort_right_float()
  rf_printer.display()
    
  print("")
  print("Unrooted average GRF:")
  grf_printer.sort_right_float()
  grf_printer.display()
    

def analyze_all():
  for datadir in os.listdir(exp.families_datasets_root):
    datadir = os.path.join(exp.families_datasets_root, datadir)
    try:
      analyze(datadir)
      print("Analyzed " + datadir)
    except:
      pass

if __name__ == '__main__':
  if (len(sys.argv) < 2):
    print("Syntax: families_dir")
    exit(1)
  print(" ".join(sys.argv))
  datadir = sys.argv[1]
  if (datadir == "all"):
    analyze_all()
  else:
    analyze(datadir)



