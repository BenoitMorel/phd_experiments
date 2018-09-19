import sys
import os

def fill_dico_and_get_total_ll(evaluator_output, dico):
  lines = open(evaluator_output).readlines()
  for line in lines[:-1]:
    split = line.split(" ")
    dico[split[0]] = float(split[3][:-1])
  return float(lines[-1].split(" ")[4][:-1])

def get_dico_sizes(dimensions_file):
  lines = open(dimensions_file).readlines()
  dico = {}
  for index in range(0, len(lines)):
    line = lines[index]
    split = line.split(" ")
    sites = int(split[0])
    taxa = int(split[1][:-1])
    dico["tree" + str(index)] = taxa * sites #[sites, taxa] 
  return dico 

if (len(sys.argv) != 5):
  print("Syntax: python compare_joint_likelihood_evaluator.py true_trees_evaluation other_trees_evaluation other_trees_name (for instance treerecs or raxml) dimensions_file")
  print("dimensions_file is a file with one line per family, witht he number of sites and the number of taxa")
  sys.exit(1)


evaluator_output_1 = sys.argv[1]
evaluator_output_2 = sys.argv[2]
second_output_name = sys.argv[3]
dimensions_file = sys.argv[4]



dico_1 = {}
dico_2 = {}

ll_1 = fill_dico_and_get_total_ll(evaluator_output_1, dico_1)
ll_2 = fill_dico_and_get_total_ll(evaluator_output_2, dico_2)

dico_sizes = get_dico_sizes(dimensions_file)

if (len(dico_1) != len(dico_2)):
  print("Error: dictionnaries have different lenghts")
  sys.exit(1)

tree_1_better = 0
tree_2_better = 0
nucleotides_tree1_better = 0
nucleotides_tree2_better = 0
nucleotides_same = 0
tree_same = 0

for key in dico_1:
  if (dico_1[key] > dico_2[key]):
    tree_1_better += 1
    nucleotides_tree1_better += dico_sizes[key]
  elif (dico_1[key] < dico_2[key]):
    tree_2_better += 1
    nucleotides_tree2_better += dico_sizes[key]
  else:
    tree_same += 1
    nucleotides_same += dico_sizes[key]

percentage_1 = float(tree_1_better) / float(len(dico_1)) * 100
percentage_2 = float(tree_2_better) / float(len(dico_1)) * 100
percentage_same = float(tree_same) / float(len(dico_1)) * 100
nucleotides_tree1_better = float(nucleotides_tree1_better) / float(tree_1_better)
nucleotides_tree2_better = float(nucleotides_tree2_better) / float(tree_2_better)
nucleotides_same = float(nucleotides_same) / float(tree_same)

print("Average size of alignments for which true trees are better: " + str(nucleotides_tree1_better))
print("Average size of alignments for which " + second_output_name + " trees are better: " + str(nucleotides_tree2_better))
print("Average size of alignments for which both are the same: " + str(nucleotides_same))
print("True trees have a better likelihood than " + second_output_name + " trees on " + str(percentage_1) + "% of the trees.")
print(second_output_name + " trees have a better likelihood than True trees on " + str(percentage_2) + "% of the trees.")
print("True trees have the same likelihood as " + second_output_name + " trees on   " + str(percentage_same) + "% of the trees.")
print("Total true     tree ll = " + str(ll_1))
print("Total " + second_output_name + " tree ll = " + str(ll_2))
if (ll_1 >= ll_2):
  print("True trees have a better total likelihood")
else:
  print("Treerecs trees have a better total likelihood")

