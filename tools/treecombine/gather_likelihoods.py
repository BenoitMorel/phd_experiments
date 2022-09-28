import os
import sys
import glob
import numpy as np

  
def extract_likelihoods(results_dir, verbose = True):
  ll1 = 0
  ll2 = 0
  improve = 0
  treated = 0
  diffs = []
  maxdiff = 0.0
  maxfam = ""
  
  for stats in glob.glob(os.path.join(results_dir, "*.stats")):
    #if ("14014_" in stats):
    #  continue
    line = open(stats).readlines()[0]
    sp = line.split()
    best_tree_number = int(float(sp[0]))
    famll1 = float(sp[1])
    famll2 = float(sp[2])
    ll1 += famll1
    ll2 += famll2
    treated += 1
    diff = famll2 - famll1
    if (diff > 0.1):
      if (verbose):
        print(os.path.basename(stats) + " " + str(diff))
      diffs.append(diff)
      improve += 1
      if (diff > maxdiff):
        maxdiff = diff
        maxfam = os.path.basename(stats)
  if (verbose):
    print("Initial ll= " + str(ll1))
    print("Final   ll= " + str(ll2))
    print("Diff      = " + str(ll2 - ll1))
    print("Av diff   = " + str(np.mean(diffs)))
    print("Median diff = " + str(np.median(diffs)))
    print("Max diff = " + str(maxdiff))
    print("Max diff = " + str(max(diffs)))
    print("Max diff for family " + maxfam)
    print("Improved families: " + str(improve) + "/" + str(treated))
  return (ll1, ll2)

if __name__ == "__main__":
  if (len(sys.argv) < 2):
    print("syntax: python " + os.path.basename(__file__) + " results_dir")
    sys.exit(1)
  results_dir = sys.argv[1]
  extract_likelihoods(results_dir)
  

