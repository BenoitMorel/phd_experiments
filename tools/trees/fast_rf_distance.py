import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
import experiments as exp
import rf_distance


def get_rf_pair(tree1, tree2):
  command = []
  command.append(exp.raxml_nompi_exec)
  command.append("--rf")
  command.append(tree1 + "," + tree2)
  
  out = subprocess.check_output(command).decode("utf-8")
  lines = out.split("\n")
  rf_abs = lines[0].split(" ")[-1]
  rf_rel = lines[1].split(" ")[-1]
  res = [float(rf_abs), float(rf_rel)]
  return res

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python fast_rf_distance.py tree1 tree2")
    sys.exit(1)
  tree1 = sys.argv[1]
  tree2 = sys.argv[2]
  


  

