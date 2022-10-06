import os
import sys
import subprocess

sys.path.insert(0, 'scripts')
import experiments as exp

def find_string_between(input_str, marker1, marker2):
  start = input_str.find(marker1) + len(marker1)
  end = input_str.find(marker2, start)
  return input_str[start:end]

def run(alignment, tree, trees, model, iqtree_prefix):
  cmd = []
  cmd.append(exp.iqtree_exec)
  cmd.append("-s")
  cmd.append(alignment)
  cmd.append("-m")
  cmd.append(model)
  cmd.append("-pre")
  cmd.append(iqtree_prefix)
  cmd.append("-redo")
  cmd.append("-z")
  cmd.append(trees)
  cmd.append("-te")
  cmd.append(tree)
  cmd.append("-n")
  cmd.append("0")
  cmd.append("-zb")
  cmd.append("10000")
  cmd.append("-zw")
  cmd.append("-au")
  cmd.append("-nt")
  cmd.append(str(40))
  #logs = subprocess.check_output(cmd, encoding='utf-8')
  print(" ".join(cmd))
  logs = subprocess.check_output(cmd)
  ll = find_string_between(logs, "BEST SCORE FOUND : ", "\n")
  return float(ll)


def run_tests(alignment, tree, trees, model, iqtree_prefix):
  run(alignment,  tree, trees, model, iqtree_prefix)




if __name__ == "__main__":
  if (len(sys.argv) < 5):
    print("syntax: python " + os.path.basename(__file__) + " alignment model tree trees iqtree_prefix]")
    sys.exit(1)
  alignment = sys.argv[1]
  tree = sys.argv[2]
  trees = sys.argv[3]
  model = sys.argv[4]
  iqtree_prefix = sys.argv[5]
  run_tests(alignment, tree, trees, model, iqtree_prefix)

  



