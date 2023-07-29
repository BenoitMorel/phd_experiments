import os
import sys
import shutil
import functools
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
sys.path.insert(0, 'tools/msa_edition')
sys.path.insert(0, 'tools/mappings')
import experiments as exp
import fam
import create_random_tree
import ete3
import read_msa
import read_tree
import get_dico
import msa_subsampler
import rename_labels

def get_string_or_none(arg):
  if (arg == "none"):
    return None
  return arg

"""
  get_families callbacks
  signature: function(alignment_dir, gene_tree_dir)
"""

def get_family_list_from_alignments(cb):
  l = []
  for f in os.listdir(cb.alignment_dir):
    if (f.startswith(".")):
      continue
    l.append(f.split(".")[0])
  return l


"""
  get_genes callbacks: return the list of genes for a given family
"""

def get_genes_from_gene_tree(cb, family):
  tree = read_tree.read_tree(cb.get_src_gene_tree(family))
  if (tree == None):
    return None
  return set(tree.get_leaf_names())

""" 
  get_species_from_gene callbacks
""" 
def get_species_from_gene_identity(gene):
  return gene

def get_species_from_gene_first_underscore(gene):
  return gene.split("_")[0]

"""
  copy_alignment callback
"""
def copy_alignment_default(cb, family): 
  src_alignment = cb.get_src_alignment(family)
  dest_alignment = cb.get_dest_alignment(family)
  shutil.copy(src_alignment, dest_alignment)

"""
  copy alignment callback, remove everything after `#` in the label names, 
  and restrict the alignment to the gene tree leaves
"""
def copy_alignment_rmcomment_prune(cb, family):
  genes_to_keep = cb.get_genes(cb, family)
  src_alignment = cb.get_src_alignment(family)
  dest_alignment = cb.get_dest_alignment(family)
  rename_labels.remove_comments(src_alignment, dest_alignment, "fasta", "#")
  msa_subsampler.prune_sequences(dest_alignment, dest_alignment, "fasta", genes_to_keep)  


"""
  main function
"""
def generate(cb):
  fam.init_top_directories(cb.datadir)

  # species tree
  if (species_tree != None):
    output_species_tree = fam.get_species_tree(cb.datadir)
    shutil.copy(cb.species_tree, output_species_tree)

  # build families one by one
  families = cb.get_family_list(cb) 
  family_number = len(families)
  family_index = 0
  for family in families:
    print(family + " (" + str(family_index) + "/" + str(family_number) + ")")
    family_index += 1
    fam.init_family_directories(cb.datadir, family)
    # gene trees
    if (cb.gene_tree_dir != None):
      src_gene_tree = cb.get_src_gene_tree(family)
      dest_gene_tree = cb.get_dest_gene_tree(family)
      shutil.copy(src_gene_tree, dest_gene_tree)
    # mappings
    genes = cb.get_genes(cb, family)
    if (genes == None):
      shutil.rmtree(fam.get_family_path(cb.datadir, family))
      print ("Family " + family + " failed!")
      continue
    # alignments
    if (alignment_dir != None):
      cb.copy_alignment(cb, family)
    species_to_gene = {}
    for gene in genes:
      species = cb.get_species_from_gene(gene)
      if (species in species_to_gene):
        species_to_gene[species].append(gene)
      else:
        species_to_gene[species] = [gene]
    mapping_file = fam.get_mappings(cb.datadir, family)
    get_dico.export_species_to_genes_to_file(species_to_gene, mapping_file)
  # final post processing
  fam.postprocess_datadir(cb.datadir)


"""
  Edit this class to customize the data processing
"""
class Callbacks():
  def __init__(self, datadir, alignment_dir, gene_tree_dir):
    self.datadir = datadir
    self.alignment_dir = alignment_dir
    self.gene_tree_dir = gene_tree_dir
    
    self.get_family_list = get_family_list_from_alignments
    self.get_genes = get_genes_from_gene_tree
    self.get_species_from_gene = get_species_from_gene_identity
    self.copy_alignment = copy_alignment_default
    self.gene_tree_suffix = ".trm.ufboot"
    self.alignment_suffix = ".trm"
    self.gene_tree_method = "true"
    self.gene_tree_model = "true"


  def get_src_gene_tree(self, family):
    return os.path.join(self.gene_tree_dir, family + self.gene_tree_suffix)

  def get_dest_gene_tree(self,family):
    return fam.get_gene_tree_path(self.datadir, family, self.gene_tree_method, self.gene_tree_model)

  def get_src_alignment(self, family):
    return os.path.join(self.alignment_dir, family + self.alignment_suffix)

  def get_dest_alignment(self, family):
    return fam.get_alignment(self.datadir, family)
  

def set_callback_hogenoms(cb):
  cb.copy_alignment = copy_alignment_rmcomment_prune
  cb.gene_tree_suffix = ".bs"
  cb.alignment_suffix = ".aln"
  callbacks.get_species_from_gene = get_species_from_gene_first_underscore  

if (__name__ == "__main__"): 
  if (len(sys.argv) < 5): 
    print("Syntax: python " + os.path.basename(__file__) + " datadir species_tree alignment_dir gene_tree_dir")
    print("To skip a parameter, use the string none")
    exit(1)
  datadir = sys.argv[1]
  species_tree = get_string_or_none(sys.argv[2])
  alignment_dir = get_string_or_none(sys.argv[3])
  gene_tree_dir = get_string_or_none(sys.argv[4])
  callbacks = Callbacks(datadir, alignment_dir, gene_tree_dir)
  
  set_callback_hogenoms(callbacks)
  
  generate(callbacks)







