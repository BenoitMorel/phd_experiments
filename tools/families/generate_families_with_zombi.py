import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "phyldog"))
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "substmodels"))
import experiments as exp
import link_file_from_gene_tree as phyldog_link
import fam
from ete3 import Tree
from ete3 import SeqGroup
from read_tree import read_tree
import model_sampler
import random

def is_ok_species_tree(output):
  species_tree = os.path.join(output, "T", "ExtantTree.nwk")
  try:
    tree = read_tree(species_tree)
    return len(tree.get_leaves() > 3)
  except:
    return False
  return False

def generate_zombi_species(species, seed, output):
  parameters_dir = os.path.join(output, "parameters")
  species_parameters_file = os.path.join(output, "SpeciesTreeParameters.tsv")
  with open(species_parameters_file, "w") as writer:
    writer.write("SPECIATION f:0.05\n")
    writer.write("EXTINCTION f:0.00\n")
    writer.write("STOPPING_RULE 1\n")
    writer.write("TOTAL_LINEAGES " + str(species) + "\n")
    writer.write("TOTAL_TIME 100\n")
    writer.write("MIN_LINEAGES 10\n")
    writer.write("MAX_LINEAGES 10000\n")
    writer.write("VERBOSE 1\n")
    writer.write("SCALE_TREE 1\n")
    writer.write("SEED " + str(seed) + "\n")
  command = []
  command.append(exp.python3())
  command.append(exp.zombi_script)
  command.append("T")
  command.append(species_parameters_file)
  command.append(output)
  print(command)
  try:
    subprocess.check_call(command)
  except:
    command[0] = exp.python3()
    subprocess.check_call(command)
  
def get_random_var_str(mean):
  return "u:0;" + str(2.0 * mean)
  #return "e:" + str(mean)

def generate_zombi_genome(families, dup_rate, loss_rate, transfer_rate, seed, output):
  parameters_dir = os.path.join(output, "parameters")
  genome_parameters_file = os.path.join(output, "GenomeTreeParameters.tsv")
  transfer_rate *= 0.01
  with open(genome_parameters_file, "w") as writer:
    writer.write("DUPLICATION " + get_random_var_str(dup_rate) + "\n")
    writer.write("LOSS " + get_random_var_str(loss_rate) + "\n")
    writer.write("TRANSFER " + get_random_var_str(transfer_rate) + "\n")
    writer.write("INVERSION f:0\n")
    writer.write("TRANSLOCATION f:0\n")
    writer.write("TRANSPOSITION f:0\n")
    writer.write("ORIGINATION f:0\n")
    writer.write("DUPLICATION_EXTENSION f:1\n")
    writer.write("TRANSFER_EXTENSION f:1\n")
    writer.write("LOSS_EXTENSION f:1\n")
    writer.write("INVERSION_EXTENSION f:0.05\n")
    writer.write("TRANSPOSITION_EXTENSION f:0.3\n")
    writer.write("TRANSLOCATION_EXTENSION f:0.3\n")
    writer.write("REPLACEMENT_TRANSFER 1\n")
    writer.write("INITIAL_GENOME_SIZE " + str(families) + "\n")
    writer.write("MIN_GENOME_SIZE 5\n")
    writer.write("GENE_LENGTH f:100\n")
    writer.write("INTERGENE_LENGTH 100\n")
    writer.write("PSEUDOGENIZATION 0.0\n")
    writer.write("SEED " + str(seed) + "\n")
    writer.write("SCALE_TREE 1\n")
    writer.write("####OUTPUT\n")
    writer.write("EVENTS_PER_BRANCH 1\n")
    writer.write("PROFILES 1\n")
    writer.write("GENE_TREES 1\n")
    writer.write("RECONCILED_TREES 0\n")
    writer.write("VERBOSE 1\n")
    writer.write("RATE_FILE False\n")
    writer.write("SCALE_RATES False\n")
    writer.write("ASSORTATIVE_TRANSFER False\n")
  command = []
  command.append(exp.python3())
  command.append(exp.zombi_script)
  command.append("Gm")
  command.append(genome_parameters_file)
  command.append(output)
  print(" ".join(command))
  subprocess.check_call(command)
  
def generate_zombi_sequence(sites, seed, output):
  parameters_dir = os.path.join(output, "parameters")
  sequence_parameters_file = os.path.join(output, "SequenceTreeParameters.tsv")
  states = ["A", "C", "G", "T"]
  with open(sequence_parameters_file, "w") as writer:
# small rate because I get very long branch lengths
    writer.write("SCALING 0.005\n") 
    
    writer.write("SEQUENCE_SIZE " + str(sites) + "\n")
    is_dna = True
    if (is_dna):
      writer.write("SEQUENCE nucleotide\n")
      tu = model_sampler.sample_from_grove("gtr", 100)
      qmatrix, symrates, freqs = tu
      for s1 in states:
        for s2 in states:
          if (s1 != s2):
            writer.write(s1 +  s2 + " " + str(qmatrix[s1 + s2]) + "\n")
      for s in states:
        writer.write(s + " " + str(freqs[s]) + "\n")
    else:
      writer.write("SEQUENCE amino-acid\n")
      writer.write("AA_MODEL LG\n")
    writer.write("ST_RATE_MULTIPLIERS n:1.0;0.2\n")
    writer.write("GF_RATE_MULTIPLIERS n:1.0;0.2\n")
    writer.write("VERBOSE 1\n")
    writer.write("SEED " + str(seed) + "\n") 
  command = []
  command.append(exp.python3())
  command.append(exp.zombi_ratecusto_script)
  command.append("S")
  command.append(sequence_parameters_file)
  command.append(output)
  subprocess.check_call(command)
  command = []
  command.append(exp.python3())
  command.append(exp.zombi_script)
  command.append("Su")
  command.append(sequence_parameters_file)
  command.append(output)
  subprocess.check_call(command)

"""
We need this renaming such that:
 - the gene name do not countain any underscore (some tools don't allow it)
 - each gene name is unique (zombi generates some genes that can have homonyms in
 other families...)
"""
def rename_with_family(zombi_gene_name, family):
  return zombi_gene_name.replace("_", "UUU") + "UUU" + family.replace("_", "UUU")


"""
  Rename all nodes from the zombi gene tree (see rename_with_family)
"""
def copy_and_rename_tree(src, dest, family):
  tree = Tree(src, 1)
  for node in tree.traverse("postorder"):
    node.name = rename_with_family(node.name, family)
  open(dest, "w").write(tree.write())


"""
  Rename all taxa in the zombi alignments (see rename_with_family)
"""
def copy_and_rename_alignment(src, dest, family):
  seqs = SeqGroup(open(src).read()) #, format="phylip_relaxed")
  new_seqs = SeqGroup() 
  for entry in seqs.get_entries():
    new_seqs.set_seq(rename_with_family(entry[0], family), entry[1])
  open(dest, "w").write(new_seqs.write())


def get_valid_families(zombi):
  families = []
  genetrees_dir = os.path.join(zombi, "G", "Gene_trees")
  for genetree_base in os.listdir(genetrees_dir):
    if (not "pruned" in genetree_base):
      continue
    genetree = os.path.join(genetrees_dir, genetree_base)
    print ("is valid " + genetree_base + "?")
    if (os.path.getsize(genetree) < 2):
      print("Not enough genes")
      continue
    family = genetree_base.split("_")[0] + "_pruned"
    alignment = os.path.join(zombi, "S", family + ".fasta")
    if (not os.path.isfile(alignment)):
      print("Invalid alignment")
      continue
    families.append(family)
  return families

"""
  Process zombi output to fit my benchmark format
"""
def zombi_to_families(zombi, datadir):
  fam.init_top_directories(datadir)
  new_ali_dir = fam.get_alignments_dir(datadir)
  new_families_dir = fam.get_families_dir(datadir)
  # species tree
  species = os.path.join(zombi, "T", "ExtantTree.nwk")
  new_species = fam.get_species_tree(datadir)
  shutil.copyfile(species, new_species)
  genetrees_dir = os.path.join(zombi, "G", "Gene_trees")
  families = get_valid_families(zombi)
  fam.init_families_directories(datadir, families)
  # species tree
  fam.convert_to_phyldog_species_tree(species, fam.get_phyldog_species_tree(datadir)) 
  for family in families:
    print("Export " + family)
    genetree = os.path.join(genetrees_dir, family + "tree.nwk")
    alignment = os.path.join(zombi, "S", family + ".fasta")
    new_family_dir = fam.get_family_path(datadir, family)
    # true trees
    print(genetree)
    copy_and_rename_tree(genetree, fam.get_true_tree(datadir, family), family)
    # alignment
    copy_and_rename_alignment(alignment, fam.get_alignment(datadir, family), family)
    # link file
    phyldog_link.generate_link_file(fam.get_true_tree(datadir, family), fam.get_mappings(datadir, family), "UUU")
  fam.postprocess_datadir(datadir)

"""
  Export gene-species mapping into decostar mapping format
  (from the per-gene phyldog mappings)
"""
def export_decostar_mappings(datadir):
  with open(fam.get_deco_mappings(datadir), "w") as writer:
    for family in fam.get_families_list(datadir):
      lines = open(fam.get_mappings(datadir, family)).readlines()
      for line in lines:
        split = line.split(":")
        species = split[0]
        split = split[1].split(";")
        for gene in split:
          writer.write(species + " " + gene.replace("\n", "") + "\n")

"""
  Export extant adjacencies from zombi output
  This will be used by decostar
  Gene names should be prefixed with species name and underscore
"""
def export_adjacencies(zombi_dir, datadir):
  species_names = Tree(fam.get_species_tree(datadir), 1).get_leaf_names()
  with open(fam.get_prefixed_adjacencies(datadir), "w") as writer:
    for species in species_names:
      genome_path = os.path.join(zombi_dir, "G", "Genomes", species + "_GENOME.tsv")
      genes = []
      for line in open(genome_path).readlines()[1:]:
        split = line.split()
        family = split[1] + "_pruned"
        name = species + "_" + split[-1]
        genes.append(rename_with_family(name, family))
      if (len(genes) <= 1):
        continue
      for i in range(0, len(genes)):
        writer.write(species + "_" + genes[i] + " " + species + "_" + genes[(i+1)%len(genes)] + "\n")
  export_decostar_mappings(datadir)

def generate_datadir(tag, species, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, seed, output):
  dirname = "zsim_" + tag  
  dirname += "_s" + str(species) + "_f" + str(families)
  dirname += "_sites" + str(sites)
  dirname += "_" + model
  dirname += "_bl" + str(bl_factor)
  dirname += "_d" + str(dup_rate) + "_l" + str(loss_rate) + "_t" + str(transfer_rate) + "_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed" + str(seed)
  datadir = os.path.join(output, dirname)
  os.makedirs(datadir)
  with open(os.path.join(datadir, "zombi_script_params.txt"), "w") as writer:
    writer.write(str(species) + " " + str(families) + " ")
    writer.write(str(sites) + " " + str(model) + " ")
    writer.write(str(bl_factor) + " " + str(dup_rate))
    writer.write(str(loss_rate) + " " + str(transfer_rate) + " " + output)
  return datadir

def generate_zombi(tag, species, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, seed, output):
  datadir = generate_datadir(tag, species, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, seed, output)
  zombi_output = os.path.join(datadir, "zombi")
  parameters_dir = os.path.join(zombi_output, "parameters")
  os.makedirs(parameters_dir)
  random.seed(seed)
  species_seed = 0
  if (not is_ok_species_tree(zombi_output)):
    generate_zombi_species(species, species_seed, zombi_output) 
    species_seed += 1000
  generate_zombi_genome(families, dup_rate, loss_rate, transfer_rate, seed, zombi_output) 
  generate_zombi_sequence(sites, seed, zombi_output)
  zombi_to_families(zombi_output, datadir)
  export_adjacencies(zombi_output, datadir)
  print("Output in: " + datadir)



if (__name__ == "__main__"): 
  if (len(sys.argv) != 12):
    print("Syntax: python3 " + os.path.basename(__file__) + " tag species families sites model bl_factor dup_rate loss_rate transfer_rate seed output")
    sys.exit(1)
  tag = sys.argv[1]
  species = int(sys.argv[2])
  families = int(sys.argv[3])
  sites = int(sys.argv[4])
  model = sys.argv[5]
  bl_factor = float(sys.argv[6])
  dup_rate = float(sys.argv[7])
  loss_rate = float(sys.argv[8])
  transfer_rate = float(sys.argv[9])
  seed = int(sys.argv[10])
  output = sys.argv[11]
  generate_zombi(tag, species, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, seed, output)
