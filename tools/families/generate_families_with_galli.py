import os
import sys
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/msa_edition')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
import ete3
import nexus_to_fasta

def remove_empty_sequences(input_msa, output_msa):
  seqs = ete3.SeqGroup(input_msa, format="fasta")
  sites = len(seqs.get_entries()[0][1])
  empty = "?" * sites
  new_seqs = ete3.SeqGroup()
  taxa = set()
  for seq in  seqs.get_entries():
    if (seq[1] != empty):
      new_seqs.set_seq(seq[0], seq[1])
      taxa.add(seq[0])
  new_seqs.write("fasta", output_msa)
  print("Number of sequences: " + str(len(new_seqs.get_entries())))
  return taxa

def generate(gene_trees_dir, ali_dir, species_tree, datadir):
  fam.init_top_directories(datadir)
  shutil.copy(species_tree, fam.get_species_tree(datadir))
  genes_to_species = {}
  for gene_tree in os.listdir(gene_trees_dir):
    family = "uce-" + gene_tree.split("-")[1]
    fam.init_family_directories(datadir, family)
    #gene_tree_src = os.path.join(gene_trees_dir, gene_tree)
    #gene_tree_list = open(gene_tree_src).readlines()
    #samples = len(gene_tree_list)
    #gene_tree_dest = fam.get_bootstrap_trees(datadir, samples, "GTR+G", family)
    #shutil.copy(gene_tree_src, gene_tree_dest)
    species_to_genes = {}
    #leaves = ete3.Tree(gene_tree_list[0]).get_leaves()
    nexus = os.path.join(ali_dir, family + ".nex")
    tempfasta = os.path.join(ali_dir, family + ".fasta")
    out_ali = fam.get_alignment(datadir, family)
    nexus_to_fasta.convert(nexus, tempfasta)
    taxa = remove_empty_sequences(tempfasta, out_ali)
    for taxon in taxa:
      species_to_genes[taxon] = [taxon]
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
  if (len(sys.argv) < 5):
    print("Syntax: python " + os.path.basename(__file__) + " gene_trees_dir ali_dir species_tree datadir")
    sys.exit(1)
  gene_trees_dir = sys.argv[1]
  ali_dir = sys.argv[2]
  species_tree = sys.argv[3]
  datadir = sys.argv[4]
  generate(gene_trees_dir, ali_dir, species_tree, datadir)

