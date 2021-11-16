import os
import sys
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
import ete3


def generate(gene_trees_dir, species_tree, datadir):
  fam.init_top_directories(datadir)
 
  genes_to_species = {}
  for gene_tree in os.listdir(gene_trees_dir):
    family = "uce-" + gene_tree.split("-")[1]
    print(family)
    fam.init_family_directories(datadir, family)
    gene_tree_src = os.path.join(gene_trees_dir, gene_tree)
    gene_tree_list = open(gene_tree_src).readlines()
    samples = len(gene_tree_list)
    gene_tree_dest = fam.get_bootstrap_trees(datadir, samples, "GTR+G", family)
    shutil.copy(gene_tree_src, gene_tree_dest)
    species_to_genes = {}
    leaves = ete3.Tree(gene_tree_list[0]).get_leaves()
    for leaf in leaves:
      species_to_genes[leaf.name] = [leaf.name]
    mapping_file = fam.get_mappings(datadir, family)
    fam.write_phyldog_mapping(species_to_genes, mapping_file)
  fam.postprocess_datadir(datadir)
      
"""
  for ali in os.listdir(ali_dir):
    family = ali.split(".")[0]
    input_ali = os.path.join(ali_dir, ali)
    output_ali = fam.get_alignment(datadir, family)
    shutil.copy(input_ali, output_ali)
    genes = get_genes(input_ali)
    species_to_genes = {}
    for gene in genes:
      species = genes_to_species[gene]
      assert (not species in species_to_genes)
      species_to_genes[species] = [gene]
    mapping_file = fam.get_mappings(datadir, family)
    fam.write_phyldog_mapping(species_to_genes, mapping_file)
  fam.postprocess_datadir(datadir)
   """ 


if (__name__ == "__main__"): 
  if (len(sys.argv) < 4):
    print("Syntax: python " + os.path.basename(__file__) + " gene_trees_dir species_tree datadir")
    sys.exit(1)
  gene_trees_dir = sys.argv[1]
  species_tree = sys.argv[2]
  datadir = sys.argv[3]
  generate(gene_trees_dir, species_tree, datadir)

