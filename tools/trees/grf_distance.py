import os
import sys
sys.path.insert(0, 'scripts')
import subprocess
import experiments as exp

def compute_distances(tree1, tree2):
  command = []
  command.append(exp.grf_eval)
  command.append(tree1)
  command.append(tree2)
  res = subprocess.check_output(command)
  lines = res.split("\n")
  rf = float(lines[0])
  grf = float(lines[1])
  print(" ".join(command))
  return (rf, grf)
  
def print_grf(tree1, tree2):
  distances = compute_distances(tree1, tree2)
  print("nRF:\t" + str(distances[0]))
  print("GRF:\t" + str(distances[1]))

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python " + os.path.basename(__file__) + " tree1 tree2")
    sys.exit(1)
  tree1 = sys.argv[1]
  tree2 = sys.argv[2]
  print_grf(tree1, tree2) 


