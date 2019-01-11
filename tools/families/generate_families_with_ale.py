import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "phyldog"))
import experiments as exp
import link_file_from_gene_tree as phyldog_link



  

def generate_ale_genome(species, families, dupRate, lossRate, transferRate, output, seed):
  generated_families = 0
  genes_dir = os.path.abspath(os.path.join(output, "gene_trees"))
  os.makedirs(genes_dir)
  oldcwd = os.getcwd()
  os.chdir(genes_dir)

  command = []
  command.append(exp.ale_simul_exec)
  command.append(str(species))
  command.append(str(species))
  command.append(str(seed))
  command.append("1")
  command.append(str(dupRate))
  command.append(str(transferRate))
  command.append(str(lossRate))
  while (len(os.listdir(genes_dir)) - 2 < families):
    subprocess.call(command) 
    print("Gene families count: " + str(len(os.listdir(genes_dir)) - 2))
  os.chdir(oldcwd)  
  gene_index = 0
  for gene in os.listdir(genes_dir):
    old = os.path.join(genes_dir, gene)
    if (gene.startswith("G")):
      new = os.path.join(genes_dir, "gene_" + str(gene_index) + ".newick")
      gene_index += 1
    if (gene.startswith("S")):
      new = os.path.join(output, "species_s.newick")
    if (gene.startswith("R")):
      new = os.path.join(output, "species_r.newick")
    shutil.move(old, new)
  
def generate_seqgen_sequence(families, sites, output, seed):
  genes_dir = os.path.abspath(os.path.join(output, "gene_trees"))
  genseq_genes_dir = os.path.abspath(os.path.join(output, "genseq_gene_trees"))
  seq_dir = os.path.abspath(os.path.join(output, "gene_sequences"))
  os.makedirs(seq_dir)
  os.makedirs(genseq_genes_dir)
  for i in range(0, len(os.listdir(genes_dir))):
    print("sequence " + str(i) + "/" + str(families))
    gene_tree = os.path.join(genes_dir, "gene_" + str(i) + ".newick")
    
    seqgen_gene_tree = os.path.join(genseq_genes_dir, "gene_" + str(i) + ".newick")
    shutil.copyfile(gene_tree, seqgen_gene_tree)
    subprocess.check_call(["sed", "-i", "s/[SDLT]//g", seqgen_gene_tree])
    sequence_file = os.path.join(seq_dir, "gene_" + str(i) + ".fasta")
    command = []
    command.append(exp.seq_gen_exec)
    command.append("-l")
    command.append(str(sites))
    command.append("-m")
    command.append("GTR")
    command.append("-of")
    command.append(seqgen_gene_tree)
    command.append("-z")
    command.append(str(int(i) + int(seed)))
    with open(sequence_file, "w") as writer:
      subprocess.check_call(command, stdout=writer)

def build_mapping_file(genetree, phyldog_mapping, treerecs_mapping):
  phyldog_writer = open(phyldog_mapping, "w")
  treerecs_writer = open(treerecs_mapping, "w")
  gene_str = open(genetree).read()
  gene_str = gene_str.replace("(", ",")
  splits = gene_str.split(",")
  dico = {}
  for split in splits:
    if (len(split) == 0):
      continue
    leaf = split.split(":")[0]
    leaf_split = leaf.split("_")
    species = leaf_split[0]
    gene = leaf
    treerecs_writer.write(gene + " " + species + "\n")
    if (not species in dico):
      dico[species] = []
    dico[species].append(gene)
  for species, genes in dico.items():
    phyldog_writer.write(species + ":" + ";".join(genes) + "\n")

def ale_to_families(ale_dir, out):
  new_ali_dir = os.path.join(out, "alignments")
  new_families_dir = os.path.join(out, "families")
  os.makedirs(new_ali_dir)
  os.makedirs(new_families_dir)
  # species tree
  species = os.path.join(ale_dir, "species_s.newick")
  new_species = os.path.join(out, "speciesTree.newick")
  shutil.copyfile(species, new_species)
  alignments_writer = open(os.path.join(out, "alignments.txt"), "w")
  old_genes_dir = os.path.join(ale_dir, "gene_trees")
  old_ali_dir = os.path.join(ale_dir, "gene_sequences")
  for genetree_base in os.listdir(old_genes_dir):
    genetree = os.path.join(old_genes_dir, genetree_base)
    if (os.path.getsize(genetree) < 2):
      continue
    family = genetree_base.split(".")[0]
    new_family = family.split("_")[1] + "_pruned"
    new_family_dir = os.path.join(new_families_dir, new_family)
    os.makedirs(new_family_dir)
    # species tree
    shutil.copyfile(new_species, os.path.join(new_family_dir, "speciesTree.newick"))
    # true trees
    shutil.copyfile(genetree, os.path.join(new_family_dir, "trueGeneTree.newick"))
    # alignment
    alignment = os.path.join(old_ali_dir, family + ".fasta")
    shutil.copyfile(alignment, os.path.join(new_family_dir, "alignment.msa"))
    shutil.copyfile(alignment, os.path.join(new_ali_dir, new_family + ".fasta"))
    alignments_writer.write(os.path.abspath(os.path.join(new_ali_dir, new_family + ".fasta")) + "\n")
    # link file
    phyldog_mapping = os.path.join(new_family_dir, "mapping.link")
    treerecs_mapping = os.path.join(new_family_dir, "treerecs_mapping.link")
    build_mapping_file(genetree, phyldog_mapping, treerecs_mapping)


def generate_ale(species, families, sites, dupRate, lossRate, transferRate, output, seed):
  dirname = "alesim_s" + str(species) + "_f" + str(families)
  dirname += "_sites" + str(sites)
  dirname += "_d" + str(dupRate) + "_l" + str(lossRate) + "_t" + str(transferRate) + "_seed" + str(seed)
  output = os.path.join(output, dirname)
  print("Writing output in " + output)
  os.makedirs(output)
  with open(os.path.join(output, "script_params.txt"), "w") as writer:
    writer.write(str(species) + " " + str(families) + " ")
    writer.write(str(sites) + " " + str(dupRate) + " ")
    writer.write(str(lossRate) + " " + str(transferRate) + " " + output)
    writer.write(" " + str(seed))
  ale_output = os.path.join(output, "ale")
  os.makedirs(ale_output)
  generate_ale_genome(species, families, dupRate, lossRate, transferRate, ale_output, seed)
  generate_seqgen_sequence(families, sites, ale_output, seed)
  print("ale output: " + ale_output)
  ale_to_families(ale_output, output)

if (len(sys.argv) != 9):
  print("Syntax: python generate_ale.py species families sites dupRate lossRate transferRate output seed")
  sys.exit(1)

species = int(sys.argv[1])
families = int(sys.argv[2])
sites = int(sys.argv[3])
dupRate = float(sys.argv[4])
lossRate = float(sys.argv[5])
transferRate = float(sys.argv[6])
output = sys.argv[7]
seed = int(sys.argv[8])


generate_ale(species, families, sites, dupRate, lossRate, transferRate, output, seed)
