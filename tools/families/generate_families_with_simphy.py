import os
import sys
import subprocess
import shutil
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "phyldog"))
sys.path.insert(0, os.path.join("tools", "trees"))
import rescale_bl
import analyze_tree
import experiments as exp

class SimphyParameters():
  def __init__(self):
    self.population = 10000
    self.speciations_per_year = 0.00001
    self.species_taxa = 30
    self.families_number = 100
    self.substitution_rate = 0.1
    self.bl = 1.0    
    self.loss_rate = 0.1
    self.dup_rate = 0.1
    self.transfer_rate = 0.0
    self.sites = 50
    self.model = "GTR"
    self.seed = 42

def build_config_file(parameters, output_dir):
  config_file = os.path.join(output_dir, "simphy_config.txt")
  with open(config_file, "w") as writer:
    writer.write("// SPECIES TREE\n")
    writer.write("-RS 1 // number of replicates\n")
    writer.write("-sb f:" + str(parameters.speciations_per_year) + "\n")
    writer.write("-sl f:" + str(parameters.species_taxa) + " // species taxa\n")
    writer.write("-su f:" + str(parameters.substitution_rate * parameters.bl * parameters.speciations_per_year) + "\n") #subsitution rate
    writer.write("-ld f:" + str(parameters.loss_rate * parameters.speciations_per_year) + "\n") # loss
    writer.write("-lb f:" + str(parameters.dup_rate * parameters.speciations_per_year) + "\n") # duplications
    writer.write("-lt f:" + str(parameters.transfer_rate * parameters.speciations_per_year) + "\n") # transfers

    writer.write("// POPULATION\n")
    writer.write("-SP f:" + str(parameters.population) + "\n")

    writer.write("// LOCUS\n")
    writer.write("-rl f:" + str(parameters.families_number) + " // locus (gene family) per replicate\n")

    writer.write("// GENERAL\n")
    writer.write("-cs " + str(parameters.seed) + "\n") 
    writer.write("-O " + str(output_dir) + " // output directory\n")
    writer.write("-OM 1 // output the mappings\n")
    writer.write("-OC 1 // log the configuration file\n")

  return config_file

def build_indelible_config_file(parameters, output_dir):
  config_file = os.path.join(output_dir, "indelible_config.txt")
  with open(config_file, "w") as writer:
    sites_mean = parameters.sites
    sites_sigma = sites_mean / 3
    writer.write("[TYPE] NUCLEOTIDE 1\n") # DNA using algorithm 1 
    writer.write("[SETTINGS] [fastaextension] fasta\n")
    writer.write("[SIMPHY-UNLINKED-MODEL] modelA \n")
    if ("GTR" == parameters.model):
      writer.write("  [submodel] GTR $(rd:2,2,2,2,2,2) // GTR with rates from a Dirichlet  \n")
      writer.write("  [statefreq] $(d:1,1,1,1)  // frequencies for T C A G sampled from a Dirichlet (1,1,1,1)\n")
    else:
      assert(False)

    #writer.write("[SIMPHY-PARTITIONS] simple [1.0 modelA $(n:" + str(sites_mean) + "," + str(sites_sigma) + ")]\n")
    writer.write("[SIMPHY-PARTITIONS] simple [1.0 modelA $(U:10," + str(int(parameters.sites) * 2) + ")]\n")

    writer.write("[SIMPHY-EVOLVE] 1 dataset \n")

  return config_file
  


def build_mapping(simphy_mapping, phyldog_mapping):
  lines = open(simphy_mapping).readlines()[1:]
  dico = {}
  for line in lines:
    if (not line.startswith("'")):
      continue
    split = line.replace("'", "").split("\t")
    gene = split[0]
    species = gene.split("_")[0]
    if (not species in dico):
      dico[species] = []
    dico[species].append(gene)
  with open(phyldog_mapping, "w") as phyldog_writer:
    for species, genes in dico.items():
      phyldog_writer.write(species + ":" + ";".join(genes) + "\n")





def run_simphy(output_dir, config_file):
  commands = []
  commands.append(exp.simphy_exec)
  commands.append("-I")
  commands.append(config_file)
  subprocess.check_call(commands)

  simphy_output_dir = os.path.join(output_dir, "1")
  species_tree = os.path.join(simphy_output_dir, "s_tree.trees")
  generations = analyze_tree.check_ultrametric_and_get_length(species_tree)
  if (False):
    rescale_bl.rescale_bl(species_tree, species_tree, 1.0 / float(generations))
    families = []
    for f in os.listdir(simphy_output_dir):
      if (f.startswith("g_trees")):
        families.append("family_" + f.split("g_trees")[1].split(".")[0])
    for family in families:
      family_number = family.split("_")[1] 
      # true trees
      gene_tree = os.path.join(simphy_output_dir, "g_trees" + family_number + ".trees")
      rescale_bl.rescale_bl(gene_tree, gene_tree, 100.0 / float(generations))
    
def run_indelible(output_dir, config_file, cores):
  commands = []
  seed = "42"
  commands.append("perl")
  commands.append(exp.simphy_indelible_wrapper)
  commands.append(output_dir)
  commands.append(config_file)
  commands.append(seed)
  commands.append(str(cores))
  subprocess.check_call(commands)

def copy_trim(input_file, output_file):
  s = open(input_file).read()
  open(output_file, "w").write(s.replace(" ", ""))

def export_to_family(output_dir):
  print("Start exporting to families format...")
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
    alignment = os.path.join(output_dir, "1", "dataset_" + family_number + ".fasta")
    
    # check that the alignment contain all characters (otherwise phyldog crashes)
    alignment_content = open(alignment).read()
    ok = True
    for c in ['A', 'C', 'G', 'T']:
      if (not c in alignment_content):
        ok = False
    if (not ok):
      shutil.rmtree(fam.get_family_path(output_dir, family))
      print("rm " + fam.get_family_path(output_dir, family))
      continue

    gene_tree = os.path.join(simphy_output_dir, "g_trees" + family_number + ".trees")
    shutil.copy(gene_tree, fam.get_true_tree(output_dir, family))
    
    simphy_mapping = os.path.join(simphy_output_dir, family_number + "l1g.maplg")
    phyldog_mapping = fam.get_mappings(output_dir, family)
    build_mapping(simphy_mapping,  phyldog_mapping)
    # alignment
    # true trees
    # alignment
    copy_trim(alignment, fam.get_alignment(output_dir, family))
    #copy_and_rename_alignment(alignment, fam.get_alignment(out, family), family)
  fam.postprocess_datadir(output_dir)

def get_output_dir(parameters, root_output):
  res = "ssim"
  res += "_s" + str(parameters.species_taxa)
  res += "_f" + str(parameters.families_number)
  res += "_sites" + str(parameters.sites)
  res += "_" + str(parameters.model).replace("+", "")
  res += "_bl" + str(parameters.bl)
  res += "_d" + str(parameters.dup_rate)
  res += "_l" + str(parameters.loss_rate)
  res += "_t" + str(parameters.transfer_rate)
  res += "_p0.0"
  res += "_pop" + str(parameters.population)
  res += "_seed" + str(parameters.seed)
  return os.path.join(root_output, res)

def generate_from_parameters(parameters, root_output):
  cores = 1
  output_dir = get_output_dir(parameters, root_output)
  exp.reset_dir(output_dir)
  config_file = build_config_file(parameters, output_dir)
  run_simphy(output_dir, config_file)
  indelible_config_file = build_indelible_config_file(parameters, output_dir)
  run_indelible(output_dir, indelible_config_file, cores)
  export_to_family(output_dir)
  print("Done! output in " + output_dir) 

def generate_simphy(species, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, perturbation, root_output, seed):
  p = SimphyParameters()
  p.species_taxa = int(species)
  p.families_number = int(families)
  p.sites = int(sites)
  p.model = model
  p.bl = float(bl_factor)
  p.dup_rate = float(dup_rate)
  p.loss_rate = float(loss_rate)
  p.transfer_rate = float(transfer_rate)
  p.seed = int(seed)
  print(perturbation)
  assert(float(perturbation) == 0.0)
  generate_from_parameters(p, root_output)


if (__name__ == "__main__"): 
  parameters = SimphyParameters()
  generate_from_parameters(parameters, exp.families_datasets_root)




