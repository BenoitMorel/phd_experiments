import os
import sys

sys.path.insert(0, 'tools/families')
import fam

def get_gene_to_species(datadir, family):
  mapping_path = fam.get_mappings(datadir, family)
  res = {}
  for line in open(mapping_path).readlines():
    split = line.replace("\n", "").split(":")
    species = split[0]
    genes = split[1].split(";")
    for gene in genes:
      res[gene] = species
  return res

def get_species_to_genes(datadir):
  res = {}
  for family in fam.get_families_list(datadir):
    mapping_path = fam.get_mappings(datadir, family)
    for line in open(mapping_path).readlines():
      split = line.replace("\n", "").split(":")
      species = split[0]
      genes = split[1].split(";")
      if (not species in res):
        res[species] = {}
      for gene in genes:
        res[species][gene] = True
  return res


