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
    translator[n[0]] = n[1]
  return translator

def get_filtered_neighbors(neighbors):
  labels = set()
  dup_labels = set()
  for n in neighbors:
    if (n[1] in labels):
      dup_labels.add(n[1])
    labels.add(n[1])
  res = set()
  for n in neighbors:
    if (not n[1] in dup_labels):
      res.add(n)
  return res


def assess_gene_tree(emf, datadir, method, subst_model, families):
  if (len(families) == 1 and families[0] == "all"):
    families = fam.get_families_list(datadir)
  per_family_hom_neighbors = find_neighbors_to_fam.get_per_family_homolog_neighbors(emf, datadir, families)
  total = 0
  for family in families:
    hom_neighbors = per_family_hom_neighbors[family]
    for neighbor_family in hom_neighbors:
      if (neighbor_family == family):
        continue
      neighbors = hom_neighbors[neighbor_family]
      neighbors = get_filtered_neighbors(neighbors)
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
      try:
        for leaf in tree1.get_leaves():
          leaf.name = translator[leaf.name]
      except:
        print("ERROR in family " + family + " and neighbor_family " + neighbor_family)
        print("Tree1 : " + tree1.write())
        print("Tree1 2 " + tree2.write())
        continue
      distance_cell = rf_distance.ete3_rf(tree1, tree2)
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

