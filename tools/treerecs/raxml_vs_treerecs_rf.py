import sys
import os

from ete3 import Tree

def print_likelihood_sums(treerecs_file):
  ale_ll = 0.0
  pll_ll = 0.0
  for line in open(treerecs_file).readlines():
    split = line.split(" ")
    for index, val in enumerate(split):
        if (val == "logLk"):
          if (split[index - 1] == "ALE"):
            ale_ll += float(split[index + 2][:-1])
          elif (split[index - 1] == "libpll"):
            pll_ll += float(split[index + 2][:-1])
  #print("Total ALE ll: " + str(ale_ll))
  #print("Total PLL ll: " + str(pll_ll))
  print("Total joint ll: " + str(pll_ll + ale_ll))

def get_relative_rf(tree1, tree2):
  rf = tree1.robinson_foulds(tree2, unrooted_trees=True)
  return float(rf[0]) / float(rf[1])

def read_list_trees(newick):
  trees = []
  for line in open(newick).readlines():
    if (line.startswith(">")):
      continue
    trees.append(Tree(line, format=1))
  return trees

def build_rf_list(trees1, trees2):
  rf_list = []
  for t1, t2 in zip(trees1, trees2):
    rf_list.append(get_relative_rf(t1, t2))
  return rf_list

def analyze_correctness(trees1, trees2, tree2_file, name):
  print("")
  print("## Analysing " + name + " trees...")
  rf_list = build_rf_list(trees1, trees2)
  rf_average = sum(rf_list) / float(len(rf_list))
  exactness_frequency = rf_list.count(0.0) / float(len(rf_list))
  print("Average relative RF with true trees: " + str(rf_average))
  print(str(exactness_frequency * 100) + "% of the trees exactly match the true trees")
  print_likelihood_sums(tree2_file)
if (len(sys.argv) != 4 and len(sys.argv) != 5):
  print("Syntax: python raxml_vs_trecs.py true_trees raxml_trees best_treerecs_trees [tree_analysis_dir]")
  sys.exit(1)

#true_trees = read_list_trees("/hits/basement/sco/morel/github/datasets/simuls/trueGeneTrees.newick")
#raxml_trees = read_list_trees("/hits/basement/sco/morel/github/datasets/simuls/geneTrees.newick")
#treerecs_trees = read_list_trees("/hits/basement/sco/morel/github/phd_experiments/results/treerecs/launch_treerecs/simuls/haswell_16/run_0/treerecs_output.newick.best")


true_trees_file = sys.argv[1]
raxml_trees_file = sys.argv[2]
treerecs_trees_file = sys.argv[3]
true_trees = read_list_trees(true_trees_file)
raxml_trees = read_list_trees(raxml_trees_file)
treerecs_trees = read_list_trees(treerecs_trees_file)
threshold_trees_dir = 0
if (len(sys.argv) == 5):
  threshold_trees_dir = sys.argv[4]

analyze_correctness(true_trees, raxml_trees, raxml_trees_file, "RAXML")
analyze_correctness(true_trees, treerecs_trees, treerecs_trees_file, "TREERECS")


if (threshold_trees_dir != 0):
  per_threshold_files = os.listdir(threshold_trees_dir)
  for per_threshold_file in per_threshold_files:
    threshold_trees = read_list_trees(os.path.join(threshold_trees_dir, per_threshold_file))
    analyze_correctness(true_trees, threshold_trees, os.path.join(threshold_trees_dir, per_threshold_file), per_threshold_file)




