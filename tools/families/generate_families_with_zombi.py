import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "phyldog"))
import experiments as exp
import link_file_from_gene_tree as phyldog_link
import fam
from ete3 import Tree
from ete3 import SeqGroup

def generate_zombi_species(species, output):
  parameters_dir = os.path.join(output, "parameters")
  species_parameters_file = os.path.join(output, "SpeciesTreeParameters.tsv")
  with open(species_parameters_file, "w") as writer:
    writer.write("SPECIATION f:0.1\n")
    writer.write("EXTINCTION f:0.02\n")
    writer.write("STOPPING_RULE 1\n")
    writer.write("TOTAL_LINEAGES " + str(species) + "\n")
    writer.write("TOTAL_TIME 1\n")
    writer.write("MIN_LINEAGES 1\n")
    writer.write("MAX_LINEAGES 10000\n")
    writer.write("VERBOSE 1\n")
    writer.write("TURNOVER F:0.0002\n")
    writer.write("LINEAGE_PROFILE 100-100;300-15000;500-50\n")
  command = []
  command.append("python3")
  command.append(exp.zombi_script)
  command.append("T")
  command.append(species_parameters_file)
  command.append(output)
  print(command)
  try:
    subprocess.check_call(command)
  except:
    command[0] = "python3.6"
    subprocess.check_call(command)
    

def generate_zombi_genome(families, dup_rate, loss_rate, transfer_rate, output):
  parameters_dir = os.path.join(output, "parameters")
  genome_parameters_file = os.path.join(output, "GenomeTreeParameters.tsv")
  with open(genome_parameters_file, "w") as writer:
    writer.write("DUPLICATION f:" + str(float(dup_rate) * float(families)) + "\n")
    writer.write("TRANSFER f:" + str(float(transfer_rate) * float(families)) + "\n")
    writer.write("LOSS f:" + str(float(loss_rate) * float(families)) + "\n")
    writer.write("INVERSION f:0\n")
    writer.write("TRANSLOCATION f:0\n")
    writer.write("ORIGINATION f:0\n")
    writer.write("DUPLICATION_EXTENSION f:1\n")
    writer.write("TRANSFER_EXTENSION f:1\n")
    writer.write("LOSS_EXTENSION f:1\n")
    writer.write("INVERSION_EXTENSION f:0.05\n")
    writer.write("TRANSLOCATION_EXTENSION f:0.3\n")
    writer.write("REPLACEMENT_TRANSFER 1\n")
    writer.write("INITIAL_GENOME_SIZE " + str(families) + "\n")
    writer.write("MIN_GENOME_SIZE 0\n")
    writer.write("GENE_LENGTH f:100\n")
    writer.write("INTERGENE_LENGTH 100\n")
    writer.write("PSEUDOGENIZATION 0.5\n")
    writer.write("####OUTPUT\n")
    writer.write("EVENTS_PER_BRANCH 1\n")
    writer.write("PROFILES 1\n")
    writer.write("GENE_TREES 1\n")
    writer.write("RECONCILED_TREES 0\n")
    writer.write("VERBOSE 1\n")

  command = []
  command.append("python3")
  command.append(exp.zombi_script)
  command.append("G")
  command.append(genome_parameters_file)
  command.append(output)
  subprocess.check_call(command)
  
def generate_zombi_sequence(sites, output):
  parameters_dir = os.path.join(output, "parameters")
  sequence_parameters_file = os.path.join(output, "SequenceTreeParameters.tsv")
  states = ["A", "C", "G", "T"]
  with open(sequence_parameters_file, "w") as writer:
# small rate because I get very long branch lengths
    writer.write("SCALING 0.005\n") 
    writer.write("SEQUENCE_SIZE " + str(sites) + "\n")
    writer.write("SEQUENCE nucleotide\n")
    for s1 in states:
      for s2 in states:
        if (s1 != s2):
          writer.write(s1 +  s2 + " 1.0\n")
    for s in states:
      writer.write(s + " 0.25\n")
    writer.write("CODON_MODEL MG\n")
    writer.write("ALPHA 1.0\n")
    writer.write("BETA 0.5\n")
    writer.write("KAPPA 1.0\n")
    writer.write("VERBOSE 1\n")
  
  command = []
  command.append("python3")
  command.append(exp.zombi_script)
  command.append("S")
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
    if (os.path.getsize(genetree) < 2):
      continue
    family = genetree_base.split("_")[0] + "_pruned"
    alignment = os.path.join(zombi, "S", family + ".fasta")
    if (not os.path.isfile(alignment)):
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
  for family in families:
    genetree = os.path.join(genetrees_dir, family + "tree.nwk")
    alignment = os.path.join(zombi, "S", family + ".fasta")
    new_family_dir = fam.get_family_path(datadir, family)
    # species tree
    fam.convert_to_phyldog_species_tree(species, fam.get_phyldog_species_tree(datadir)) 
    # true trees
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

def generate_datadir(species, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, output):
  dirname = "zsim_s" + str(species) + "_f" + str(families)
  dirname += "_sites" + str(sites)
  dirname += "_" + model
  dirname += "_bl" + str(bl_factor)
  dirname += "_d" + str(dup_rate) + "_l" + str(loss_rate) + "_t" + str(transfer_rate)
  datadir = os.path.join(output, dirname)
  os.makedirs(datadir)
  with open(os.path.join(datadir, "zombi_script_params.txt"), "w") as writer:
    writer.write(str(species) + " " + str(families) + " ")
    writer.write(str(sites) + " " + str(model) + " ")
    writer.write(str(bl_factor) + " " + str(dup_rate))
    writer.write(str(loss_rate) + " " + str(transfer_rate) + " " + output)
  return datadir

def generate_zombi(species, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, output):
  datadir = generate_datadir(species, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, output)
  zombi_output = os.path.join(datadir, "zombi")
  parameters_dir = os.path.join(zombi_output, "parameters")
  os.makedirs(parameters_dir)
  generate_zombi_species(species, zombi_output) 
  generate_zombi_genome(families, dup_rate, loss_rate, transfer_rate, zombi_output) 
  generate_zombi_sequence(sites, zombi_output)
  zombi_to_families(zombi_output, datadir)
  export_adjacencies(zombi_output, datadir)
  print("Output in: " + datadir)



if (__name__ == "__main__"): 
  if (len(sys.argv) != 10):
    print("Syntax: python generate_zombi.py species families sites model bl_factor dup_rate loss_rate transfer_rate output")
    sys.exit(1)
  species = int(sys.argv[1])
  families = int(sys.argv[2])
  sites = int(sys.argv[3])
  model = sys.argv[4]
  bl_factor = float(sys.argv[5])
  dup_rate = float(sys.argv[6])
  loss_rate = float(sys.argv[7])
  transfer_rate = float(sys.argv[8])
  output = sys.argv[9]
  generate_zombi(species, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, output)
