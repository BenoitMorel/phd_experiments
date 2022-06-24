import os
import sys
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
sys.path.insert(0, 'tools/mappings')
import read_synteny_from_emf as emf_reader
import get_dico
import read_tree
import fam

def get_neighbor_gene_name(gene_name, info, after):
  if (not gene_name in info.allgenes):
    return None
  emf_gene = info.allgenes[gene_name]
  neighbor = None
  if (after):
    neighbor = emf_gene.after
  else:
    neighbor = emf_gene.before
  return neighbor

def get_neighbor_family(gene_name, info, gene_to_family, after):
  neighbor = get_neighbor_gene_name(gene_name, info, after)
  if (neighbor == None or not neighbor in gene_to_family):
    return "NONE"
  
  
  res1 = gene_to_family[neighbor]
  res2 = "None"
  
  double_neighbor = get_neighbor_gene_name(neighbor, info, after)
  if (double_neighbor == None or not double_neighbor in gene_to_family):
    res2 = "None"
  else:
    res2 = gene_to_family[double_neighbor]

  if (after):
    return res1 + "-" + res2
  else:
    return res2 + "-" + res1
  

def print_gene_tree(emf, datadir, family, method, subst_model):
  print("Reading emf...")
  info = emf_reader.read(emf)
  print("Reading all mappings")
  gene_to_family = get_dico.get_gene_to_family(datadir)
  gene_names = get_dico.get_genes(datadir, family)
  gene_tree_path = fam.build_gene_tree_path(datadir, subst_model, family, method)
  gene_tree = read_tree.read_tree(gene_tree_path)
  gene_to_species = get_dico.get_gene_to_species(datadir, family)
  
  for leaf in gene_tree.get_leaves():
    species = gene_to_species[leaf.name]
    family_before = get_neighbor_family(leaf.name, info, gene_to_family, False)
    family_after = get_neighbor_family(leaf.name, info, gene_to_family, True)
    leaf.name = family_before + "-" + species + "-" + family_after
    
  print(gene_tree.write())


if (__name__ == "__main__"):
  if (len(sys.argv) < 6):
    print("Syntax python " + os.path.basename(__file__) + " emf datadir family method subst_model")
  emf = sys.argv[1]
  datadir = sys.argv[2]
  family = sys.argv[3]
  method = sys.argv[4]
  subst_model = sys.argv[5]
  print_gene_tree(emf, datadir, family, method, subst_model)

