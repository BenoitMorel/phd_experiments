import sys
import os
from ete3 import Tree


def get_rf(tree1, tree2):
  return tree1.robinson_foulds(tree2, unrooted_trees=True)[0]

def get_relative_rf(tree1, tree2):
  rf = tree1.robinson_foulds(tree2, unrooted_trees=True)
  return float(rf[0]) / float(rf[1])



def analyse(dataset_dir, pargenes_dir):
  analysed_msas = 0
  total_raxml_rrf = 0.0
  total_treerecs_rrf = 0.0
  total_jointsearch_rrf = 0.0
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
    raxml_rrf = get_relative_rf(true_tree, raxml_tree)
    treerecs_rrf = get_relative_rf(true_tree, treerecs_tree)
    jointsearch_rrf = get_relative_rf(true_tree, jointsearch_tree)

    best_rrf = min(raxml_rrf, treerecs_rrf, jointsearch_rrf)
    if (best_rrf == raxml_rrf):
      best_raxml += 1
    if (best_rrf == treerecs_rrf):
      best_treerecs += 1
    if (best_rrf == jointsearch_rrf):
      best_jointsearch += 1

    #print(msa + " \t " + str(raxml_rf) + " ")
    total_raxml_rrf += raxml_rrf
    total_treerecs_rrf += treerecs_rrf
    total_jointsearch_rrf += jointsearch_rrf
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
  print("Number of gene families for which a method reaches the smallest relative RF to the true trees compared with the other methods:")
  print("- Raxml:       " + str(best_raxml))
  print("- Treerecs:    " + str(best_treerecs))
  print("- JointSearch: " + str(best_jointsearch))

if __name__ == '__main__':
  if (len(sys.argv) != 3):
    print("Syntax: families_dir pargenes_dir")
    exit(1)
  dataset_dir = sys.argv[1]
  pargenes_dir = sys.argv[1]
  analyse(dataset_dir, pargenes_dir)



