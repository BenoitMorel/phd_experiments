import math
import os
import sys
import subprocess
import shutil
import fam
import copy
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "phyldog"))
sys.path.insert(0, os.path.join("tools", "trees"))
import rescale_bl
import analyze_tree
import experiments as exp
import discordance_rate
import sample_missing_data

class SimphyParameters():
  def __init__(self):
    self.tag = "tag"
    self.prefix = "ssim"
    self.population = 10000
    self.speciations_per_year = 0.000000005
    self.species_taxa = 25
    self.families_number = 200
    self.bl = 1.0    
    self.loss_rate = 1.0
    self.dup_rate = 1.0
    self.transfer_rate = 0.0
    self.sites = 50
    self.model = "GTR"
    self.seed = 42
    self.distance_hgt = False

def build_config_file(parameters, output_dir):
  config_file = os.path.join(output_dir, "simphy_config.txt")
  with open(config_file, "w") as writer:
    writer.write("// SPECIES TREE\n")
    # number of replicates
    writer.write("-RS 1 // number of replicates\n")
    # speciation rates (speciation per yer)
    writer.write("-sb f:" + str(parameters.speciations_per_year) + "\n")
    # number of species taxa 
    writer.write("-sl f:" + str(parameters.species_taxa) + "\n")
    # species tree height in years (I don't understand this)
    writer.write("-st ln:21.25,0.2\n")
    # substitution rate
    writer.write("-su ln:-21.9," + str(0.1 * parameters.bl) + "\n")
    # L, D, T global rates 
    lognormal_scale = 1.0
    lognormal_location = 0.0 #math.log(1.0 - 0.5 * pow(lognormal_scale, 2.0))
    lognormal_mean = math.exp((lognormal_location + pow(lognormal_scale, 2.0)) / 2.0)
    loss_freq =     0.00000000049 * parameters.loss_rate / lognormal_mean
    dup_freq =      0.00000000049 * parameters.dup_rate / lognormal_mean
    transfer_freq = 0.00000000049 * parameters.transfer_rate / lognormal_mean
    assert(loss_freq == dup_freq)
    writer.write("-gd f:" + str(loss_freq) + "\n")
    writer.write("-gb f:" + str(dup_freq) + "\n")
    writer.write("-gt f:" + str(transfer_freq) +"\n")

    # L, D, T per family rates
    writer.write("-ld sl:" + str(lognormal_location) + "," + str(lognormal_scale) + ",gd\n")
    writer.write("-lb f:ld\n")
    writer.write("-lt sl:" + str(lognormal_location) + "," + str(lognormal_scale) + ",gt\n")
    #writer.write("-lt f:ld\n")
    
    lk = 0
    if (parameters.distance_hgt):
        lk = 1
    writer.write("-lk " + str(lk) + "\n")
    
    #writer.write("-lb sl:" + str(lognormal_location) + "," + str(lognormal_scale) + ",gb\n")
    #writer.write("-lt sl:" + str(lognormal_location) + "," + str(lognormal_scale) + ",gt\n")

    writer.write("// POPULATION\n")
    writer.write("-SP f:" + str(parameters.population) + "\n")

    writer.write("// LOCUS\n")
    writer.write("-rl f:" + str(parameters.families_number) + " // locus (gene family) per replicate\n")

    writer.write("// Subsitution rates heterogeneity parameters\n")
    writer.write("-hs ln:1.5,1\n")
    writer.write("-hl ln:1.551533,0.6931472\n")
    writer.write("-hg ln:1.5,1\n")


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
    sites_min = str(20)
    sites_max = str(2 * int(parameters.sites) - 20)
    #sites_sigma = sites_mean / 3
    writer.write("[TYPE] NUCLEOTIDE 1\n") # DNA using algorithm 1 
    writer.write("[SETTINGS] [fastaextension] fasta\n")
    writer.write("[SIMPHY-UNLINKED-MODEL] modelA \n")
    if ("GTR" == parameters.model):
      writer.write("  [submodel] GTR $(rd:16,3,5,5,6,15) // GTR with rates from a Dirichlet  \n")
      writer.write("  [statefreq] $(d:36,26,28,32)  // frequencies for T C A G sampled from a Dirichlet \n")
      writer.write("[rates] 0 $(e:2) 0 // Site-specific rate heterogeneities: 0 p-inv, alpha from an E(2) and using a continuous gamma distribution.\n")
    else:
      assert(False)

    #writer.write("[SIMPHY-PARTITIONS] simple [1.0 modelA $(f:" + str(int(parameters.sites)) + ")]\n")
    writer.write("[SIMPHY-PARTITIONS] simple [1.0 modelA $(u::" + sites_min + "," + sites_max + ")]\n")

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

def export_to_family(output_dir, replicate = 1):
  print("Start exporting to families format...")
  fam.init_top_directories(output_dir)
  simphy_output_dir = os.path.join(output_dir, str(replicate))
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
  res = parameters.prefix
  res += "_" + parameters.tag
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
  res += "_mu" + str(parameters.mu)
  res += "_theta" + str(parameters.theta)
  """
  res += "_hgt" 
  if (parameters.distance_hgt):
    res += "dist"
  else:
    res += "unif"
  """
  res += "_seed" + str(parameters.seed)
  return os.path.join(root_output, res)

def compute_and_write_discordance_rate(parameters, output_dir):
  d = 0.0
  if (int(parameters.population) > 20):
      if (parameters.dup_rate == 0.0 and parameters.transfer_rate == 0.0):
        d = discordance_rate.get_discordance_rate(output_dir)
      else:
        no_dtl_parameters = copy.deepcopy(parameters) 
        no_dtl_parameters.dup_rate = 0.0
        no_dtl_parameters.loss_rate = 0.0
        no_dtl_parameters.transfer_rate = 0.0
        temp_output_dir = generate_from_parameters(no_dtl_parameters, output_dir)
        d = fam.get_discordance_rate(temp_output_dir)
        shutil.rmtree(temp_output_dir)
    
  print("Discordance rate: " + str(d))
  fam.write_discordance_rate(output_dir, d)


def generate_from_parameters(parameters, root_output):
  cores = 1
  output_dir = get_output_dir(parameters, root_output)
  os.mkdir(output_dir)
  exp.reset_dir(output_dir)
  config_file = build_config_file(parameters, output_dir)
  run_simphy(output_dir, config_file)
  indelible_config_file = build_indelible_config_file(parameters, output_dir)
  run_indelible(output_dir, indelible_config_file, cores)
  export_to_family(output_dir)
  compute_and_write_discordance_rate(parameters, output_dir)
  print("Done! output in " + output_dir) 
  return output_dir

def generate_simphy(tag, species, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, perturbation, population, mu, theta, root_output, seed):
  p = SimphyParameters()
  p.tag = tag
  p.species_taxa = int(species)
  p.families_number = int(families)
  p.sites = int(sites)
  p.model = model
  p.bl = float(bl_factor)
  p.dup_rate = float(dup_rate)
  p.loss_rate = float(loss_rate)
  p.transfer_rate = float(transfer_rate)
  p.seed = int(seed)
  p.population = population
  p.mu = float(mu)
  p.theta = float(theta)
  print(perturbation)
  assert(float(perturbation) == 0.0)
  if (p.mu != 1.0 or p.theta != 0.0):
    p.prefix = p.prefix + "temp"
    temp_output_dir = generate_from_parameters(p, root_output)
    p.prefix = p.prefix[:-4]
    output_dir = get_output_dir(p, root_output) 
    print("Now move " + temp_output_dir + " to " + output_dir)
    sample_missing_data.sample_missing_data(temp_output_dir, output_dir, p.mu, p.theta)
    shutil.rmtree(temp_output_dir)
  else:
    generate_from_parameters(p, root_output)
   

if (__name__ == "__main__"): 
  parameters = SimphyParameters()
  generate_from_parameters(parameters, exp.families_datasets_root)




