import os
import sys
import subprocess
import shutil
import hashlib
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "phyldog"))
sys.path.insert(0, os.path.join("tools", "trees"))
import experiments as exp
import link_file_from_gene_tree as phyldog_link
import sequence_model
import rescale_bl
import analyze_tree
from ete3 import Tree
from ete3 import SeqGroup


def generate_jprime_species(species, output, seed):
  print("species tree...")
  species_parameters_file = os.path.join(output, "SpeciesTreeParameters.tsv")
  command = []
  command.append("java")
  command.append("-jar")
  command.append(exp.jprime_jar)
  command.append("HostTreeGen")
  command.append("-nox")
  command.append(str(species))
  command.append("-s")
  command.append(str(seed))
  command.append("-min")
  command.append("5")
  command.append("1")
  command.append("0")
  command.append(os.path.join(output, "species"))
  subprocess.check_call(command)
  
  species_tree = os.path.join(output, "species.pruned.tree")
  subprocess.check_call(["sed", "-i", "s/1.443047701658439E-4/0.0001443047701658439/g", species_tree])
  print(open(os.path.join(output, "species.pruned.tree")).read()) 

def generate_jprime_genome(families, dup_rate, loss_rate, transfer_rate, output, seed):
  for i in range(0, families):
    print("gene " + str(i) + "/" + str(families))
    command = []
    command.append("java")
    command.append("-jar")
    command.append(exp.jprime_jar)
    command.append("GuestTreeGen")
    #command.append("-nox")
    command.append("-max")
    command.append("5000000")
    command.append("-s")
    command.append(str(i + seed))
    command.append(os.path.join(output, "species.pruned.tree"))
    command.append(str(dup_rate))
    command.append(str(loss_rate))
    command.append(str(transfer_rate))
    command.append(os.path.join(output, str(i) + "_gene"))
    subprocess.check_call(command)
    

  
def generate_seqgen_sequence(families, sites, model, output, seed):
  model_samples = sequence_model.get_model_samples(model)
  seqgen_model_cmd_samples = [sequence_model.model_to_seqgen_cmd(m) for m in model_samples]
  for i in range(0, families):
    print("sequence " + str(i) + "/" + str(families))
    jprime_tree = os.path.join(output, str(i) + "_gene.pruned.tree")
    sequence_file = os.path.join(output, str(i) + ".fasta")
    command = []
    command.append(exp.seq_gen_exec)
    command.append("-l")
    command.append(str(sites))
    command.append("-of")
    command.append(jprime_tree)
    command.append("-z")
    command.append(str(int(i) + int(seed)))
    command.extend(seqgen_model_cmd_samples[i % len(seqgen_model_cmd_samples)])
    with open(sequence_file, "w") as writer:
      subprocess.check_call(command, stdout=writer)

def build_mapping_file(jprime_mapping, phyldog_mapping, treerecs_mapping, family):
  phyldog_writer = open(phyldog_mapping, "w")
  treerecs_writer = open(treerecs_mapping, "w")
  lines = open(jprime_mapping).readlines()
  dico = {}
  for line in lines:
    split = line.split("\t")
    split[1] = split[1].rstrip()
    split[0] = rename_with_family(split[0], family)
    treerecs_writer.write(split[0] + " " + split[1] + "\n")
    if (not split[1] in dico):
      dico[split[1]] = []
    dico[split[1]].append(split[0])
  for species, genes in dico.items():
    phyldog_writer.write(species + ":" + ";".join(genes) + "\n")

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


def jprime_to_families(jprime, out):
  fam.init_top_directories(out)
  families = []
  for genetree_base in os.listdir(jprime):
    if (not "gene.pruned.tree" in genetree_base):
      continue
    genetree = os.path.join(jprime, genetree_base)
    if (os.path.getsize(genetree) < 2):
      continue
    family_number = genetree_base.split("_")[0]
    family =  family_number + "_pruned"
    families.append(family)
  fam.init_families_directories(out, families)
  # species tree
  species = os.path.join(jprime, "species.pruned.tree")
  shutil.copyfile(species, fam.get_species_tree(out))
  for family in families:
    family_number = family.split("_")[0] 
    # true trees
    gene_tree = os.path.join(jprime, family_number + "_gene.pruned.tree")
    # alignment
    alignment = os.path.join(jprime, family_number + ".fasta")
    # true trees
    copy_and_rename_tree(gene_tree, fam.get_true_tree(out, family), family)
    # alignment
    copy_and_rename_alignment(alignment, fam.get_alignment(out, family), family)
    # link file
    jprime_mapping = os.path.join(jprime, family_number + "_gene.pruned.leafmap")
    phyldog_mapping = fam.get_mappings(out, family)
    treerecs_mapping = fam.get_treerecs_mappings(out, family)
    build_mapping_file(jprime_mapping, phyldog_mapping, treerecs_mapping, family)

def rescale_trees(jprime_output, families, bl_factor):
  for i in range(0, families):
    tree = os.path.join(jprime_output, str(i) + "_gene.pruned.tree")
    subprocess.check_call(["sed", "-i", "s/\[[^][]*\]//g", tree])
    subprocess.check_call(["sed", "-i", "s/)[^:]*:/):/g", tree])
    rescale_bl.rescale_bl(tree, tree, bl_factor)
    subprocess.check_call(["sed", "-i", "s/)1:/):/g", tree])

def get_output(species, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, perturbation):
  dirname = "jsim"
  if (float(transfer_rate) != 0.0):
    dirname += "dtl"
  dirname += "_s" + str(species) + "_f" + str(families)
  dirname += "_sites" + str(sites)
  dirname += "_" + model
  dirname += "_bl" + str(bl_factor)
  dirname += "_d" + str(dup_rate) + "_l" + str(loss_rate) 
  dirname += "_t" + str(transfer_rate)  
  dirname += "_p" + str(perturbation)  
  
  return dirname

def generate_jprime(species, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, perturbation, root_output, seed):
  to_hash = str(species) + str(families) + str(sites) + model + str(bl_factor) + str(dup_rate) + str(loss_rate) + str(transfer_rate) + str(seed) + str(perturbation)
  md5 = hashlib.md5(to_hash.encode())
  output = os.path.join(root_output, "jprime_temp_" + str(md5.hexdigest()))
  shutil.rmtree(output, True)
  print("Writing output in " + output)
  os.makedirs(output)
  jprime_output = os.path.join(output, "jprime")
  os.makedirs(jprime_output)
  with open(os.path.join(jprime_output, "jprime_script_params.txt"), "w") as writer:
    writer.write(str(species) + " " + str(families) + " ")
    writer.write(str(sites) + " " + str(model) + " ")
    writer.write(str(bl_factor)+ " " + str(dup_rate) + " ")
    writer.write(str(loss_rate) + " " + str(transfer_rate) + " " + str(perturbation) + " " + output)
    writer.write(" " + str(seed))
  generate_jprime_species(species, jprime_output, seed) 
  generate_jprime_genome(families, dup_rate, loss_rate, transfer_rate, jprime_output, seed)
  rescale_trees(jprime_output, families, bl_factor)
  generate_seqgen_sequence(families, sites, model, jprime_output, seed)
  print("jprime output: " + jprime_output)
  jprime_to_families(jprime_output, output)
  
  species_nodes = analyze_tree.get_tree_taxa_number(os.path.join(jprime_output, "species.pruned.tree"))
  new_output = os.path.join(root_output, get_output(species_nodes, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, perturbation))
  shutil.move(output, new_output)
  fam.postprocess_datadir(new_output)
  fam.perturbate_species_tree(new_output, perturbation)
  
  print("Final output directory: " + new_output)
  print("")

if (__name__ == "__main__"): 

  if (len(sys.argv) != 12 or not (sys.argv[4] in sequence_model.get_model_sample_names())):
    if (len(sys.argv) != 12):
      print("Invalid number of parameters")
    print("Syntax: python generate_jprime.py species_time_interval families sites model bl_scaler dup_rate loss_rate transfer_rate species_perturbation output seed")
    print("model should be one of " + str(sequence_model.get_model_sample_names()))
    exit(1)

  species = int(sys.argv[1])
  families = int(sys.argv[2])
  sites = int(sys.argv[3])
  model = sys.argv[4]
  bl_factor = float(sys.argv[5])
  dup_rate = float(sys.argv[6])
  loss_rate = float(sys.argv[7])
  transfer_rate = float(sys.argv[8])
  perturbation = float(sys.argv[9])
  output = sys.argv[10]
  seed = int(sys.argv[11])

  generate_jprime(species, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, perturbation, output, seed)
