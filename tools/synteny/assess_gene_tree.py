import os
import sys
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import find_neighbors_to_fam
import fam
import read_tree
import rf_distance
import prune

def get_induced_gene_tree(datadir, family, method, subst_model, leaf_set):
  gene_tree_path = fam.get_gene_tree(datadir, subst_model, family, method)
  gene_tree = read_tree.read_tree(gene_tree_path)
  #gene_tree.prune(leaf_set)
  return prune.fast_prune(gene_tree, leaf_set)
  #return gene_tree

def get_translator(neighbors):
  translator = {}
  for n in neighbors:
    translator[n[1]] = n[0]
  return translator

def assess_gene_tree(emf, datadir, method, subst_model, families):
  if (len(families) == 1 and families[0] == "all"):
    families = fam.get_families_list(datadir)
  all_neighbors = find_neighbors_to_fam.find_neighbors(emf, datadir, families)
  total = 0
  for family in families:
    for direction in [0, 1]:
      neighbor_families = all_neighbors[family][direction]
      for neighbor_family in neighbor_families:
        if (neighbor_family == family):
          continue
        neighbors = neighbor_families[neighbor_family]
        if (len(neighbors) < 4):
          continue
        
        leaf_set1 = set()
        leaf_set2 = set()
        for neighbor in neighbors:
          leaf_set1.add(neighbor[0])
          leaf_set2.add(neighbor[1])
        tree1 = get_induced_gene_tree(datadir, family, method, subst_model, leaf_set1)
        tree2 = get_induced_gene_tree(datadir, neighbor_family, method, subst_model, leaf_set2)
        translator = get_translator(neighbors)
        for leaf in tree2.get_leaves():
          leaf.name = translator[leaf.name]
        distance_cell = rf_distance.ete3_rf(tree1, tree2)
        print("absolute rf = " + str(distance_cell[0]))
        total += distance_cell[0]
  print("Total: " + str(total))




if (__name__ == "__main__"):
  if (len(sys.argv) < 6):
    print("Syntax python " + os.path.basename(__file__) + " emf datadir method subst_model families")
  emf = sys.argv[1]
  datadir = sys.argv[2]
  method = sys.argv[3]
  subst_model = sys.argv[4]
  families = sys.argv[5:]
  assess_gene_tree(emf, datadir, method, subst_model, families)

