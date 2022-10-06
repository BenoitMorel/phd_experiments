import os
import sys
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/mappings')
sys.path.insert(0, 'tools/msa_edition')
import fam
import get_dico
from read_msa import read_msa
    



def get_leaf_names(ali):    
  names = set()
  for entry in ali.get_entries():   
    names.add(entry[0])
  return names
        

def clean_mappings(datadir):
  for family in fam.get_families_list(datadir):
    ali = read_msa(fam.get_alignment(datadir, family))
    genes = get_leaf_names(ali) 
    old_dico = get_dico.get_gene_to_species(datadir, family)
    new_dico = {}
    for gene in genes:
      species = old_dico[gene]
      if (not species in new_dico):
        new_dico[species] = []
      new_dico[species].append(gene)
    phyldog = fam.get_mappings(datadir, family)
    treerecs = fam.get_treerecs_mappings(datadir, family)
    get_dico.export_species_to_genes_to_file(new_dico, phyldog)
    fam.convert_phyldog_to_treerecs_mapping(phyldog, treerecs)


if (__name__== "__main__"):
  if len(sys.argv) != 2:
    print("Syntax: python " + os.path.basename(__file__) + " datadir ")
    sys.exit(1)
  datadir = sys.argv[1]
  clean_mappings(datadir)

