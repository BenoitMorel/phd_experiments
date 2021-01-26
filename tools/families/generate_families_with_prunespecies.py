import sys
import os
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
sys.path.insert(0, 'tools/mappings')
import experiments as exp
import fam
from ete3 import Tree
from ete3 import SeqGroup
import get_dico



def get_species_to_keep(inputdir, species_to_prune):
  species_tree = Tree(fam.get_species_tree(inputdir), format = 1)
  all_species = set(species_tree.get_leaf_names())
  return all_species - set(species_to_prune)


def extract_species_tree(inputdir, outputdir, species_to_keep):
  tree = Tree(fam.get_species_tree(inputdir), format = 1)
  tree.prune(species_to_keep)
  tree.write(outfile = fam.get_species_tree(outputdir))

def extract_family(family, inputdir, outputdir, gene_method, subst_model, species_to_keep):
  tree = Tree(fam.get_true_tree(inputdir, family), format = 1)
  species_to_genes = get_dico.get_species_to_genes_family(inputdir, family)
  new_species_to_genes = {}
  genes_to_keep = set()
  for species in species_to_genes:
    if (species in species_to_keep):
      new_species_to_genes[species] = species_to_genes[species]
      genes_to_keep |= set(species_to_genes[species])
  if (len(genes_to_keep) < 4):
    return
  fam.init_family_directories(outputdir, family)
  fam.write_phyldog_mapping(new_species_to_genes, fam.get_mappings(outputdir, family))
  input_gene_tree = fam.build_gene_tree_path(inputdir, subst_model, family, gene_method)
  if (os.path.isfile(input_gene_tree)):
    tree = Tree(input_gene_tree, format = 1)
    tree.prune(genes_to_keep)
    output_gene_tree = fam.build_gene_tree_path(outputdir, subst_model, family, gene_method)
    tree.write(outfile = output_gene_tree)
  input_ali = fam.get_alignment(inputdir, family)
  if (os.path.isfile(input_ali)):
    msa = SeqGroup(input_ali)
    new_msa = SeqGroup()
    for entry in msa:
      if (entry[0] in genes_to_keep):
        new_msa.set_seq(entry[0], entry[1])
    new_msa.write(outfile = fam.get_alignment(outputdir, family))

def generate(inputdir, outputdir, gene_method, subst_model, species_to_prune):
  print("Species to prune: " + " ".join(species_to_prune))
  fam.init_top_directories(outputdir)   
  with open(os.path.join(fam.get_misc_dir(outputdir), "info.txt"), "w") as writer:
    writer.write("Extracted from " + os.path.basename(inputdir))
    writer.write(" by removing the species:\n" + "\n".join(species_to_prune))
  species_to_keep = get_species_to_keep(inputdir, species_to_prune)
  print("Species to keep: " + " ".join(species_to_keep))
  extract_species_tree(inputdir, outputdir, species_to_keep)
  families = fam.get_families_list(inputdir)
  for family in families:
    extract_family(family, inputdir, outputdir, gene_method, subst_model, species_to_keep)
  fam.postprocess_datadir(outputdir)
  print("Result datadir in " + outputdir)

  

if (__name__ == "__main__"): 
  if (len(sys.argv) < 5): 
    print("Syntax: python " + os.path.basename(__file__) + " input output gene_method subst_model species1 [species2 species3 ...]")
    exit(1)
  inputdir = sys.argv[1]
  outputdir = sys.argv[2]
  gene_method = sys.argv[3]
  subst_model = sys.argv[4]
  species_to_prune = sys.argv[5:]
  generate(inputdir, outputdir, gene_method, subst_model, species_to_prune)



