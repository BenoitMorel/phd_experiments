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
  print(input_ali)
  msa = None
  try:
    msa	= SeqGroup(input_ali, format = "phylip_relaxed")
  except:
    msa	= SeqGroup(input_ali, format = "iphylip_relaxed")
  genes = set()
  for entry in msa.get_entries():
    genes.add(entry[0])
  return genes
  
def generate(msa_dir, genetree_dir, species_tree, outputdir):
  fam.init_top_directories(outputdir)
  output_species_tree = fam.get_species_tree(outputdir)
  shutil.copy(species_tree, output_species_tree)
  all_species = set(get_species(species_tree)) 

  for ali in os.listdir(msa_dir):
    family = ali.replace(".phy", "")
    print(family)
    fam.init_family_directories(outputdir, family)
    input_ali = os.path.join(msa_dir, ali)
    output_ali = fam.get_alignment(outputdir, family)
    shutil.copy(input_ali, output_ali)
    genes = get_genes(input_ali)
    species_to_genes = {}
    for gene in genes:
      species = gene.split("_")[0]
      if ("merge" in family):
        species = gene.split("_")[1]
      if (not species in all_species):
        print(family)
        print("error with gene " + gene + ": species " + species  + " is not in the species tree")
      assert (species in all_species)
      if (not species in species_to_genes):
        species_to_genes[species] =  []
      species_to_genes[species].append(gene)
    mapping_file = fam.get_mappings(outputdir, family)
    fam.write_phyldog_mapping(species_to_genes, mapping_file)
    input_gene_tree = os.path.join(genetree_dir, family + ".treelist")
    output_gene_tree = fam.build_gene_tree_path(outputdir, "GTR", family, "phylobayes") 
    shutil.copy(input_gene_tree, output_gene_tree)
  fam.postprocess_datadir(outputdir)
    


if (__name__ == "__main__"): 
  if (len(sys.argv) < 5):
    print("Syntax: python" + os.path.basename(__file__) + " msa_dir genetree_dir species_tree output_dir")
    sys.exit(1)
  msa_dir = sys.argv[1]
  genetree_dir = sys.argv[2]
  species_tree = sys.argv[3]
  outputdir = sys.argv[4]
  generate(msa_dir, genetree_dir, species_tree, outputdir)

