import os
import sys
import subprocess
import shutil
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "phyldog"))
sys.path.insert(0, os.path.join("tools", "trees"))
import experiments as exp

def build_config_file(output_dir):
  species_taxa = 10
  families_number = 5
  config_file = os.path.join(output_dir, "simphy_config.txt")
  with open(config_file, "w") as writer:
    writer.write("// SPECIES TREE\n")
    writer.write("-RS 1 // number of replicates\n")
    writer.write("-sb f:0.00001 // speciations per year\n")
    writer.write("-sl f:" + str(species_taxa) + " // species taxa\n")

    writer.write("-lb f:0.000005\n") # duplications

    writer.write("// POPULATION\n")
    writer.write("-SP f:1000 // population size\n")

    writer.write("// LOCUS\n")
    writer.write("-rl f:" + str(families_number) + " // locus (gene family) per replicate\n")

    writer.write("// GENERAL\n")
    writer.write("-cs 42 // random seed\n")
    writer.write("-O " + str(output_dir) + " // output directory\n")
    writer.write("-OM 1 // output the mappings\n")
    writer.write("-OC 1 // log the configuration file\n")

  return config_file

def build_mapping(simphy_mapping, phyldog_mapping):
  lines = open(simphy_mapping).readlines()[1:]
  dico = {}
  for line in lines:
    split = line.replace("'", "").split("\t")
    gene = split[0]
    species = gene.split("_")[0]
    if (not species in dico):
      dico[species] = []
    dico[species].append(gene)
  with open(phyldog_mapping, "w") as phyldog_writer:
    for species, genes in dico.items():
      phyldog_writer.write(species + ":" + ";".join(genes) + "\n")

def run_simphy(config_file):
  commands = []
  commands.append(exp.simphy_exec)
  commands.append("-I")
  commands.append(config_file)
  subprocess.check_call(commands)

def export_to_family(output_dir):
  fam.init_top_directories(output_dir)
  simphy_output_dir = os.path.join(output_dir, "1")
  families = []
  for f in os.listdir(simphy_output_dir):
    if (f.startswith("g_trees")):
      families.append("family_" + f.split("g_trees")[1].split(".")[0])
  fam.init_families_directories(output_dir, families)
  # species tree
  species = os.path.join(simphy_output_dir, "s_tree.trees")
  shutil.copyfile(species, fam.get_species_tree(output_dir))
  for family in families:
    family_number = family.split("_")[1] 
    # true trees
    gene_tree = os.path.join(simphy_output_dir, "g_trees" + family_number + ".trees")
    shutil.copy(gene_tree, fam.get_true_tree(output_dir, family))
    
    simphy_mapping = os.path.join(simphy_output_dir, family_number + "l1g.maplg")
    phyldog_mapping = fam.get_mappings(output_dir, family)
    build_mapping(simphy_mapping,  phyldog_mapping)
    # alignment
    #alignment = os.path.join(, family_number + ".fasta")
    # true trees
    # alignment
    #copy_and_rename_alignment(alignment, fam.get_alignment(out, family), family)
  

output_dir = "../BenoitDatasets/families/simphy_test"
exp.reset_dir(output_dir)
config_file = build_config_file(output_dir)
run_simphy(config_file)
export_to_family(output_dir)






