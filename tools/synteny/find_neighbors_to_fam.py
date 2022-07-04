import os
import sys
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/mappings')
import read_synteny_from_emf as emf_reader
import get_dico




"""
  return all the neighbors of the genes in the given family,
  as a set of tuples (gene name, its neighbor name)
"""
def get_neighbors(info, gene_to_family, datadir, family, before = True, after = True):
  gene_names = get_dico.get_genes(datadir, family)
  neighbors = set()
  for name in gene_names:
    if (not name in info.allgenes):
      continue
    emf_gene = info.allgenes[name]
    if (after):
      neighbor = emf_gene.after
      if (neighbor != None and neighbor in gene_to_family):
        neighbors.add((emf_gene.gene, neighbor))
    if (before):
      neighbor = emf_gene.before
      if (neighbor != None and neighbor in gene_to_family):
        neighbors.add((emf_gene.gene, neighbor))
  return neighbors


"""
  returns homolog_neighbors such that homolog_neighbors[f]
  is the set of tuples (gene, neighbor_gene) such that gene belongs
  to family "family" and neighbor_gene is a neighbor of gene and
  belongs to family "f"
"""
def get_homolog_neighbors(info, gene_to_family, datadir, family, before = True, after = True):
  neighbors = get_neighbors(info, gene_to_family, datadir, family, before, after)
  homolog_neighbors = {}
  
  for neighbor in neighbors:
    neighbor_fam = gene_to_family[neighbor[1]]
    if (not neighbor_fam in homolog_neighbors):
      homolog_neighbors[neighbor_fam] = []
    homolog_neighbors[gene_to_family[neighbor[1]]].append(neighbor)
  return homolog_neighbors

def get_per_family_homolog_neighbors(emf, datadir, families, before = True, after = True):
  print("Reading emf...")
  info = emf_reader.read(emf)
  print("Reading gene to family dict")
  gene_to_family = get_dico.get_gene_to_family(datadir)
  print("Filling the per-family homolog neighbors...")
  res = {}
  for family in families:
    res[family] = get_homolog_neighbors(info, gene_to_family, datadir, family, before, after)  
  return res


def find_neighbors_fam(info, gene_to_family, datadir, family, after = True, verbose = False):
  gene_names = get_dico.get_genes(datadir, family)
  print("Number of genes in the family " + family + " " + str(len(gene_names)))
  neighbors = {}
  for name in gene_names:
    if (not name in info.allgenes):
      continue
    emf_gene = info.allgenes[name]
    if (emf_gene == None):
      print("Gene " + name + " is not present in the emf file")
      continue
    neighbor = emf_gene.after
    if (not after):
      neighbor = emf_gene.before
    if (neighbor == None):
      continue
    if (not neighbor in gene_to_family):
      continue
    family_neighbor = gene_to_family[neighbor]
    
    if (not family_neighbor in neighbors):
      neighbors[family_neighbor] = []
    neighbors[family_neighbor].append((emf_gene.gene, neighbor))

  if (verbose):
    for family in neighbors:
      print(family + ": " +  str(len(neighbors[family])))
    print("")
    print("")
  return neighbors






def find_neighbors(emf, datadir, families, verbose = False):
  print("Reading emf...")
  info = emf_reader.read(emf)
  gene_to_family = get_dico.get_gene_to_family(datadir)
  neighbors = {}
  for family in families:
    n1 = find_neighbors_fam(info, gene_to_family, datadir, family, True, verbose)
    n2 = find_neighbors_fam(info, gene_to_family, datadir, family, False, verbose)
    neighbors[family] = [n1, n2]
  return neighbors

if (__name__ == "__main__"):
  if (len(sys.argv) < 4):
    print("Syntax python " + os.path.basename(__file__) + " emf datadir family")
  emf = sys.argv[1]
  datadir = sys.argv[2]
  families = sys.argv[3:]
  find_neighbors(emf, datadir, families, verbose = True)
