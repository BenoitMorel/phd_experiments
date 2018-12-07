import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "phyldog"))
import experiments as exp
import link_file_from_gene_tree as phyldog_link



def generate_jprime_species(species, output):
  print("species tree...")
  species_parameters_file = os.path.join(output, "SpeciesTreeParameters.tsv")
  command = []
  command.append("java")
  command.append("-jar")
  command.append(exp.jprime_jar)
  command.append("HostTreeGen")
  command.append("-nox")
  command.append(str(species))
  command.append("1")
  command.append("0")
  command.append(os.path.join(output, "species"))
  subprocess.check_call(command)
  print(open(os.path.join(output, "species.pruned.tree")).read()) 

def generate_jprime_genome(families, dupRate, lossRate, transferRate, output):
  for i in range(0, families):
    print("gene " + str(i) + "/" + str(families))
    command = []
    command.append("java")
    command.append("-jar")
    command.append(exp.jprime_jar)
    command.append("GuestTreeGen")
    #command.append("-nox")
    command.append("-max")
    command.append("500")
    command.append("-s")
    command.append(str(i))
    command.append(os.path.join(output, "species.pruned.tree"))
    command.append(str(dupRate))
    command.append(str(lossRate))
    command.append(str(transferRate))
    command.append(os.path.join(output, str(i) + "_gene"))
    subprocess.check_call(command)
    

  
def generate_seqgen_sequence(families, sites, output):
  for i in range(0, families):
    print("sequence " + str(i) + "/" + str(families))

    jprime_tree = os.path.join(output, str(i) + "_gene.pruned.tree")
    seqgene_tree = os.path.join(output, str(i) + "_gene.seqgen.tree")
    sequence_file = os.path.join(output, str(i) + ".fasta")
    shutil.copyfile(jprime_tree, seqgene_tree)
    subprocess.check_call(["sed", "-i", "s/\[[^][]*\]//g", jprime_tree])
    subprocess.check_call(["sed", "-i", "s/)[^:]*:/):/g", jprime_tree])
    subprocess.check_call(["sed", "-i", "s/\[[^][]*\]//g", seqgene_tree])
    subprocess.check_call(["sed", "-i", "s/)[^:]*:/):/g", seqgene_tree])
    command = []
    command.append(exp.seq_gen_exec)
    command.append("-l")
    command.append(str(sites))
    command.append("-m")
    command.append("GTR")
    command.append("-of")
    command.append(seqgene_tree)
    with open(sequence_file, "w") as writer:
      subprocess.check_call(command, stdout=writer)

def build_mapping_file(gprime_mapping, phyldog_mapping, treerecs_mapping):
  phyldog_writer = open(phyldog_mapping, "w")
  treerecs_writer = open(treerecs_mapping, "w")
  lines = open(gprime_mapping).readlines()
  dico = {}
  for line in lines:
    split = line.split("\t")
    treerecs_writer.write(split[0] + " " + split[1])
    split[1] = split[1][:-1]
    if (not split[1] in dico):
      dico[split[1]] = []
    dico[split[1]].append(split[0])
  for species, genes in dico.items():
    phyldog_writer.write(species + ":" + ";".join(genes) + "\n")

def gprime_to_families(gprime, out):
  new_ali_dir = os.path.join(out, "alignments")
  new_families_dir = os.path.join(out, "families")
  os.makedirs(new_ali_dir)
  os.makedirs(new_families_dir)
  # species tree
  species = os.path.join(gprime, "species.pruned.tree")
  new_species = os.path.join(out, "speciesTree.newick")
  shutil.copyfile(species, new_species)
  alignments_writer = open(os.path.join(out, "alignments.txt"), "w")
  for genetree_base in os.listdir(gprime):
    if (not "gene.pruned.tree" in genetree_base):
      continue
    genetree = os.path.join(gprime, genetree_base)
    if (os.path.getsize(genetree) < 2):
      continue
    family_number = genetree_base.split("_")[0]
    family =  family_number + "_pruned"
    print(family)
    new_family_dir = os.path.join(new_families_dir, family)
    os.makedirs(new_family_dir)
    # species tree
    shutil.copyfile(new_species, os.path.join(new_family_dir, "speciesTree.newick"))
    # true trees
    shutil.copyfile(genetree, os.path.join(new_family_dir, "trueGeneTree.newick"))
    # alignment
    alignment_base = family_number + ".fasta" 
    alignment = os.path.join(gprime, alignment_base)
    shutil.copyfile(alignment, os.path.join(new_family_dir, "alignment.msa"))
    shutil.copyfile(alignment, os.path.join(new_ali_dir, alignment_base))
    alignments_writer.write(os.path.abspath(os.path.join(new_ali_dir, alignment_base)) + "\n")
    # link file
    gprime_mapping = os.path.join(gprime, genetree_base[:-4] + "leafmap")
    phyldog_mapping = os.path.join(new_family_dir, "mapping.link")
    treerecs_mapping = os.path.join(new_family_dir, "treerecs_mapping.link")
    build_mapping_file(gprime_mapping, phyldog_mapping, treerecs_mapping)
    #phyldog_link.generate_link_file(genetree, os.path.join(new_family_dir, "mapping.link"), "_")


def generate_jprime(species, families, sites, dupRate, lossRate, transferRate, output):
  dirname = "jsim_s" + str(species) + "_f" + str(families)
  dirname += "_sites" + str(sites)
  dirname += "_d" + str(dupRate) + "_l" + str(lossRate)
  output = os.path.join(output, dirname)
  print("Writing output in " + output)
  os.makedirs(output)
  with open(os.path.join(output, "jprime_script_params.txt"), "w") as writer:
    writer.write(str(species) + " " + str(families) + " ")
    writer.write(str(sites) + " " + str(dupRate) + " ")
    writer.write(str(lossRate) + " " + str(transferRate) + " " + output)
  jprime_output = os.path.join(output, "jprime")
  os.makedirs(jprime_output)
  generate_jprime_species(species, jprime_output) 
  generate_jprime_genome(families, dupRate, lossRate, transferRate, jprime_output)
  generate_seqgen_sequence(families, sites, jprime_output)
  print("jprime output: " + jprime_output)
  gprime_to_families(jprime_output, output)

if (len(sys.argv) != 8):
  print("Syntax: python generate_jprime.py species_time_interval families sites dupRate lossRate transferRate output")
  sys.exit(1)

species = int(sys.argv[1])
families = int(sys.argv[2])
sites = int(sys.argv[3])
dupRate = float(sys.argv[4])
lossRate = float(sys.argv[5])
transferRate = float(sys.argv[6])
output = sys.argv[7]

generate_jprime(species, families, sites, dupRate, lossRate, transferRate, output)
