import os
import sys
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/mappings')
import read_synteny_from_emf as emf_reader
import get_dico




"""
  return all the neighbors (as gene names) of the genes in the given family
"""
def get_neighbors(info, gene_to_family, datadir, family, before = True, after = True):
  gene_names = get_dico.get_genes(datadir, family)
  neighbors = []
  for name in gene_names:
    if (not name in info.allgenes):
      continue
    emf_gene = info.allgenes[name]
    if (after):
      neighbor = emf_gene.after
      if (neighbor != None and neighbor in gene_to_family):
        neighbors.append(gene_to_family[neighbor])
  return neighbors



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
