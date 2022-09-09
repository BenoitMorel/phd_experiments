import os
import sys
import glob
import numpy as np

def extract_likelihoods(results_dir, hard_only = True):
  ll1 = 0
  ll2 = 0
  improve = 0
  treated = 0
  diffs = []
  for stats in glob.glob(os.path.join(results_dir, "*.stats")):
    line = open(stats).readlines()[0]
    sp = line.split()
    best_tree_number = int(float(sp[0]))
    if (not hard_only or best_tree_number < 2):
      diff1 = float(sp[1])
      diff2 = float(sp[2])
      ll1 += diff1
      ll2 += diff2
      treated += 1
      if (diff2 - diff1 > 0.000001):
        diffs.append(diff2 - diff1)
        #print(str(diff2 - diff1) + " " + str(diff1))
        improve += 1
  print("Initial ll= " + str(ll1))
  print("Final   ll= " + str(ll2))
  print("Diff      = " + str(ll2 - ll1))
  print("Av diff   = " + str(np.mean(diffs)))
  print("Median diff = " + str(np.median(diffs)))
  print("Max diff = " + str(np.max(diffs)))
  print("Improved families: " + str(improve) + "/" + str(treated))

if __name__ == "__main__":
  if (len(sys.argv) < 2):
    print("syntax: python " + os.path.basename(__file__) + " results_dir")
    sys.exit(1)
  results_dir = sys.argv[1]
  extract_likelihoods(results_dir, False)
  #extract_likelihoods(results_dir, True)
  

