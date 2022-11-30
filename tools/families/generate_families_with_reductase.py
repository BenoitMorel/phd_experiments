import os
import sys
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
import read_tree
from ete3 import SeqGroup


def get_species(species_tree):
  tree = read_tree.read_tree(species_tree)
  return tree.get_leaf_names()

def get_genes(input_ali):
  msa = SeqGroup(input_ali)
  genes = set()
  for entry in msa.get_entries():
    genes.add(entry[0])
  return genes
  
def generate(inputdir, outputdir):
  fam.init_top_directories(outputdir)
  ali_dir = os.path.join(inputdir, "bmge")
  input_species_tree = os.path.join(inputdir, "S.tree")
  output_species_tree = fam.get_species_tree(outputdir)
  shutil.copy(input_species_tree, output_species_tree)
  all_species = set(get_species(input_species_tree)) 

  for ali in os.listdir(ali_dir):
    family = ali.split(".")[0]
    fam.init_family_directories(outputdir, family)
    input_ali = os.path.join(ali_dir, ali)
    output_ali = fam.get_alignment(outputdir, family)
    shutil.copy(input_ali, output_ali)
    genes = get_genes(input_ali)
    species_to_genes = {}
    for gene in genes:
      species = gene.split("_")[0]
      assert (species in all_species)
      if (not species in species_to_genes):
        species_to_genes[species] =  []
      species_to_genes[species].append(gene)
    mapping_file = fam.get_mappings(outputdir, family)
    fam.write_phyldog_mapping(species_to_genes, mapping_file)
  fam.postprocess_datadir(outputdir)
    


if (__name__ == "__main__"): 
  if (len(sys.argv) < 3):
    print("Syntax: python" + os.path.basename(__file__) + " input_dir output_dir")
  inputdir = sys.argv[1]
  outputdir = sys.argv[2]
  generate(inputdir, outputdir)

