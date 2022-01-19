import sys
import os
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
sys.path.insert(0, 'tools/mappings')
sys.path.insert(0, 'tools/msa_edition')
import experiments as exp
import fam
from ete3 import Tree
from ete3 import SeqGroup
import get_dico
import read_msa
import ete3

import re


def is_ok(seq, max_gap_ratio, is_dna):
  pattern = "[ACGTacgt]"
  if (not is_dna):
    pattern = "[ACDEFGHIKLMNPQRSTVWYabcdefghiklmnpqrstvwy]"
  l = float(len(seq))
  c = float(len(re.findall(pattern, seq)))
  return c/l >= (1.0 - max_gap_ratio)
  
def get_filtered_msa(msa, max_gap_ratio, is_dna):
  new_msa = SeqGroup()
  for entry in msa.get_entries():
    if (is_ok(entry[1], max_gap_ratio, is_dna)):
      gene = entry[0]
      if (gene == "Smilax_bona"):
        gene = "Smilax_bona-nox"
      if (gene == "Riccia_sp"):
        gene = "Ricciocarpos_natans"
      new_msa.set_seq(gene, entry[1])
  return new_msa

def generate(inputdir, max_gap_ratio, is_dna):
  outputdir = os.path.normpath(inputdir) +  "_maxgapratio" + str(max_gap_ratio)
  print(outputdir)
  fam.init_top_directories(outputdir)   
  ok_species = set()
  total_genes_before = 0
  total_genes_after = 0
  for family in fam.get_families_list(inputdir):
    msa = read_msa.read_msa(fam.get_alignment(inputdir, family))
    total_genes_before += len(msa.get_entries())
    new_msa = get_filtered_msa(msa, max_gap_ratio, is_dna)
    if (len(new_msa.get_entries()) < 4):
      continue
    total_genes_after += len(new_msa.get_entries())
    fam.init_family_directories(outputdir, family)
    new_msa.write(format = "fasta", outfile = fam.get_alignment(outputdir, family)) 
    species_to_genes = {}
    old_gene_to_species = get_dico.get_gene_to_species(inputdir, family)
    for entry in new_msa.get_entries():
      gene = entry[0]
      species = old_gene_to_species[gene]
      ok_species.add(species)
      if (not species in species_to_genes):
        species_to_genes[species] = []
      species_to_genes[species].append(gene)
    fam.write_phyldog_mapping(species_to_genes, fam.get_mappings(outputdir, family))
  print("Initial families: " + str(len(fam.get_families_list(inputdir))))
  print("Initial genes: " + str(total_genes_before))
  print("Final   families: " + str(len(fam.get_families_list(outputdir))))
  print("Final   genes: " + str(total_genes_after))
  species_tree = ete3.Tree(fam.get_species_tree(inputdir), format = 1)
  species_tree.prune(ok_species)
  species_tree.write(format = 1, outfile = fam.get_species_tree(outputdir))
  fam.postprocess_datadir(outputdir)



if (__name__ == "__main__"): 
  if (len(sys.argv) < 3): 
    print("Syntax: python " + os.path.basename(__file__) + " datadir max_gap_ratio is_dna (1 or 0)")
    exit(1)
  datadir = sys.argv[1]
  max_gap_ratio = float(sys.argv[2])
  is_dna = int(sys.argv[3]) == 1
  generate(datadir, max_gap_ratio, is_dna)


