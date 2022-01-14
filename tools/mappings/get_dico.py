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

def get_genes(datadir, family):
  mapping_path = fam.get_mappings(datadir, family)
  res = set()
  for line in open(mapping_path).readlines():
    split = line.replace("\n", "").split(":")
    genes = split[1].split(";")
    for gene in genes:
      res.add(gene)
  return res


def export_species_to_genes_to_file(species_to_genes, output):
  with open(output, "w") as writer: 
    for species in species_to_genes:
      genes = list(species_to_genes[species])
      if (len(genes) > 0):
          writer.write(species)
          writer.write(":")
          writer.write(";".join(genes))
          writer.write("\n")

def export_species_to_genes(species_to_genes, datadir, family):
  export_species_to_genes_to_file(species_to_genes, fam.get_mappings(datadir, family))


def get_species_to_genes(datadir):
  res = {}
  for family in fam.get_families_list(datadir):
    mapping_path = fam.get_mappings(datadir, family)
    for line in open(mapping_path).readlines():
      split = line.replace("\n", "").split(":")
      species = split[0]
      genes = split[1].split(";")
      if (not species in res):
        #res[species] = {}
        res[species] = set()
      for gene in genes:
        #res[species][gene] = True
        res[species].add(gene)
  return res

def get_species_to_genes_family(datadir, family):
  res = {}
  mapping_path = fam.get_mappings(datadir, family)
  for line in open(mapping_path).readlines():
    split = line.replace("\n", "").split(":")
    species = split[0]
    genes = split[1].split(";")
    if (not species in res):
      res[species] = []
    for gene in genes:
      res[species].append(gene)
  return res


