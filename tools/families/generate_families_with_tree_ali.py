"""
  All input families should be single-copy families
  ali_dir: gene alignments named family.msa
  gene_tree_dir: gene trees named family.newick
"""


import os
import sys
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
from ete3 import Tree



def fill_dict_from_gene_tree(datadir, family, gene_tree_path):
  leaves = Tree(gene_tree_path).get_leaves()
  species_to_genes = {}
  for leaf in leaves:
    species_to_genes[leaf.name] = [leaf.name]
  mapping_file = fam.get_mappings(datadir, family)
  fam.write_phyldog_mapping(species_to_genes, mapping_file)



def generate(ali_dir, gene_tree_dir, species_tree, datadir):
  fam.init_top_directories(datadir)
  shutil.copyfile(species_tree, fam.get_species_tree(datadir))
  for ali_name in os.listdir(ali_dir):
    family = ali_name.split(".")[0]
    ali = os.path.join(ali_dir, ali_name)
    gene_tree = os.path.join(gene_tree_dir, family + ".newick")
    if (not os.path.isfile(gene_tree)):
      print("Skipping family " + family)
      continue
    fam.init_family_directories(datadir, family)
    shutil.copyfile(ali, fam.get_alignment(datadir, family))
    shutil.copyfile(gene_tree, fam.get_true_tree(datadir, family))
    fill_dict_from_gene_tree(datadir, family, gene_tree)

  fam.postprocess_datadir(datadir)

if (__name__ == "__main__"): 
  if (len(sys.argv) < 5):
    print("Syntax: python" + os.path.basename(__file__) + " ali_dir gene_tree_dir species_tree outputdir")
    sys.exit(1)
  ali_dir = sys.argv[1]
  gene_tree_dir = sys.argv[2]
  species_tree = sys.argv[3]
  datadir = sys.argv[4]
  generate(ali_dir, gene_tree_dir, species_tree, datadir)


