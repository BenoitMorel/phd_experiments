import sys
import os

def fill_dico_and_get_total_ll(evaluator_output, dico):
  lines = open(evaluator_output).readlines()
  for line in lines[:-1]:
    split = line.split(" ")
    dico[split[0]] = float(split[3][:-1])
  return float(lines[-1].split(" ")[4][:-1])


if (len(sys.argv) != 3):
  print("Syntax: python compare_joint_likelihood_evaluator.py true_trees_evaluation treerecs_trees_evaluation")
  sys.exit(1)


evaluator_output_1 = sys.argv[1]
evaluator_output_2 = sys.argv[2]

dico_1 = {}
dico_2 = {}

ll_1 = fill_dico_and_get_total_ll(evaluator_output_1, dico_1)
ll_2 = fill_dico_and_get_total_ll(evaluator_output_2, dico_2)

if (len(dico_1) != len(dico_2)):
  print("Error: dictionnaries have different lenghts")
  sys.exit(1)

tree_1_better = 0
tree_2_better = 0
tree_same = 0

for key in dico_1:
  if (dico_1[key] > dico_2[key]):
    tree_1_better += 1
  elif (dico_1[key] < dico_2[key]):
    tree_2_better += 1
  else:
    tree_same += 1

percentage_1 = float(tree_1_better) / float(len(dico_1)) * 100
percentage_2 = float(tree_2_better) / float(len(dico_1)) * 100
percentage_same = float(tree_same) / float(len(dico_1)) * 100

print("True trees have a better likelihood than Treerecs trees on " + str(percentage_1) + "% of the trees.")
print("Treerecs trees have a better likelihood than True trees on " + str(percentage_2) + "% of the trees.")
print("True trees have the same likelihood as Treerecs trees on   " + str(percentage_same) + "% of the trees.")
print("ll_1: " + str(ll_1))
print("ll_2: " + str(ll_2))
if (ll_1 >= ll_2):
  print("True trees have a better total likelihood")
else:
  print("Treerecs trees have a better total likelihood")

