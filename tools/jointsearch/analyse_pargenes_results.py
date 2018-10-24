import sys
import os
from ete3 import Tree


def get_rf(tree1, tree2):
  return tree1.robinson_foulds(tree2, unrooted_trees=True)[0]

def get_relative_rf(tree1, tree2):
  rf = tree1.robinson_foulds(tree2, unrooted_trees=True)
  return float(rf[0]) / float(rf[1])

def get_nodes(tree1):
  rf = tree1.robinson_foulds(tree1, unrooted_trees=True)
  return rf[1]



def analyse(dataset_dir, pargenes_dir):
  analysed_msas = 0
  total_raxml_rrf = 0.0
  total_treerecs_rrf = 0.0
  total_jointsearch_rrf = 0.0
  total_raxml_rf = 0.0
  total_treerecs_rf = 0.0
  total_jointsearch_rf = 0.0
  total_nodes_number = 0
  true_raxml = 0
  true_treerecs = 0
  true_jointsearch = 0
  best_raxml = 0
  best_treerecs = 0
  best_jointsearch = 0
  for msa in os.listdir(dataset_dir):   
    family_path = os.path.join(dataset_dir, msa)
    try:
      true_tree = Tree(os.path.join(family_path, "trueGeneTree.newick"), format=1) 
      raxml_tree = Tree(os.path.join(family_path, "raxmlGeneTree.newick"), format=1)
      treerecs_tree = Tree(os.path.join(family_path, "treerecsGeneTree.newick"), format=1)
      jointsearch_tree = Tree(os.path.join(pargenes_dir, "results", msa, "jointsearch.newick"), format=1)
    
    except:
      continue
    #raxml_rf = get_rf(true_tree, raxml_tree)
    raxml_rf_cell = raxml_tree.robinson_foulds(true_tree, unrooted_trees=True)
    treerecs_rf_cell = treerecs_tree.robinson_foulds(true_tree, unrooted_trees=True)
    jointsearch_rf_cell = jointsearch_tree.robinson_foulds(true_tree, unrooted_trees=True)
    
    raxml_rrf = float(raxml_rf_cell[0]) / float(raxml_rf_cell[1])
    treerecs_rrf = float(treerecs_rf_cell[0]) / float(treerecs_rf_cell[1])
    jointsearch_rrf = float(jointsearch_rf_cell[0]) / float(jointsearch_rf_cell[1])

    raxml_rf = float(raxml_rf_cell[0])
    treerecs_rf = float(treerecs_rf_cell[0])
    jointsearch_rf = float(jointsearch_rf_cell[0])
    
    total_nodes_number += raxml_rf_cell[1]

    best_rrf = min(raxml_rrf, treerecs_rrf, jointsearch_rrf)
    if (best_rrf == raxml_rrf):
      best_raxml += 1
    if (best_rrf == treerecs_rrf):
      best_treerecs += 1
    if (best_rrf == jointsearch_rrf):
      best_jointsearch += 1

    total_raxml_rrf += raxml_rrf
    total_treerecs_rrf += treerecs_rrf
    total_jointsearch_rrf += jointsearch_rrf
    total_raxml_rf += raxml_rf
    total_treerecs_rf += treerecs_rf
    total_jointsearch_rf += jointsearch_rf

    if (raxml_rf_cell[0] == 0):
      true_raxml += 1
    if (treerecs_rf_cell[0] == 0):
      true_treerecs += 1
    if (jointsearch_rf_cell[0] == 0):
      true_jointsearch += 1
    analysed_msas += 1

  if (analysed_msas == 0):
    print("did not manage to analyse any MSA")
    exit(1)
  print("Number of gene families: " + str(analysed_msas))
  print("")
  print("Average (over the gene families) relative RF distance to the true trees:")
  print("- Raxml-ng:   " + str(total_raxml_rrf / float(analysed_msas)))
  print("- Treerecs:   " + str(total_treerecs_rrf / float(analysed_msas)))
  print("- JointSearch " + str(total_jointsearch_rrf / float(analysed_msas)))
  print("")
  print("Normalized average (over the gene families) RF distance to the true trees:")
  print("- Raxml-ng:   " + str(total_raxml_rf / float(total_nodes_number)))
  print("- Treerecs:   " + str(total_treerecs_rf / float(total_nodes_number)))
  print("- JointSearch " + str(total_jointsearch_rf / float(total_nodes_number)))
  print("")
  print("Number of gene families for which a method reaches the smallest relative RF to the true trees compared with the other methods:")
  print("- Raxml:       " + str(best_raxml))
  print("- Treerecs:    " + str(best_treerecs))
  print("- JointSearch: " + str(best_jointsearch))
  print("")
  print("Number of gene families for which a method finds the true tree:")
  print("- Raxml:       " + str(true_raxml))
  print("- Treerecs:    " + str(true_treerecs))
  print("- JointSearch: " + str(true_jointsearch))

if __name__ == '__main__':
  if (len(sys.argv) != 3):
    print("Syntax: families_dir pargenes_dir")
    exit(1)
  dataset_dir = sys.argv[1]
  pargenes_dir = sys.argv[2]
  analyse(dataset_dir, pargenes_dir)



