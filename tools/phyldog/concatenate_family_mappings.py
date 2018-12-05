import os
import sys

def concatenate_family_mappings(families_dir, output_file):
  species_genes = {}  
  for family in os.listdir(families_dir):
    lines = []
    try:
      lines = open(os.path.join(families_dir, family, "phyldog", "phyldogMapping.link")).readlines()
    except:
      continue
    for line in lines:
      split = line.split(":")
      species = split[0]
      genes_to_add = split[1][:-1].split(";")
      if not species in species_genes:
        species_genes[species] = []
      genes = species_genes[species]
      for gene in genes_to_add:
        genes.append(gene)
  with open(output_file, "w") as writer:
    for species, genes in species_genes.iteritems():
      writer.write(species + ":" + ",".join(genes) + "\n")
    

  

if (len(sys.argv) != 3):
  print("Syntax: python concatenate_family_mappings.py families_dir output_file" )
  exit(1)

families_dir = sys.argv[1]
output_file_str = sys.argv[2]
concatenate_family_mappings(families_dir, output_file_str)
