import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import perturbate

def mkdir(directory):
  try:
    os.makedirs(directory)
  except:
    pass

######################
# global directories
######################

def get_datasets_family_path():
  return os.path.join(exp.benoit_datasets_root, "families")

def get_datadir(dataset):
  return os.path.abspath(os.path.join(exp.benoit_datasets_root, "families", dataset))

def get_species_dir(datadir):
  return os.path.join(datadir, "species_trees")

def get_adjacencies_dir(datadir):
  return os.path.join(datadir, "adjacencies")

def get_families_dir(datadir):
  return os.path.join(datadir, "families")

def get_misc_dir(datadir):
  return os.path.join(datadir, "misc")

def get_alignments_dir(datadir):
  return os.path.join(datadir, "alignments")

def get_run_dir(datadir, subst_model = None, run_name = None):
  res =  os.path.join(datadir, "runs")
  if (subst_model != None):
    res = os.path.join(res, subst_model)
  if (run_name != None):
    res = os.path.join(res, run_name)
  return res

######################
# global files
#####################

def get_species_tree(datadir, subst_model= None, method = "true"):
  if (method == "random"):
    return "random"
  elif (method == "NJ"):
    return "NJ"
  if (method == "true"):
    return os.path.join(get_species_dir(datadir), "speciesTree.newick")
  if (subst_model != None):
    return os.path.join(get_species_dir(datadir), method + "." + subst_model + ".speciesTree.newick")
  else:
    return os.path.join(get_species_dir(datadir), method + ".speciesTree.newick")

def get_true_species_tree(datadir):
  return os.path.join(get_species_dir(datadir), "trueSpeciesTree.newick")

def get_phyldog_species_tree(datadir):
  return os.path.join(get_species_dir(datadir), "phyldogSpeciesTree.newick")

def get_adjacencies(datadir):
  return os.path.join(get_adjacencies_dir(datadir), "adjacencies.txt")

def get_prefixed_adjacencies(datadir):
  return get_adjacencies(datadir) + ".prefixed"

def get_deco_mappings(datadir):
  return os.path.join(get_adjacencies_dir(datadir), "deco_mappings.txt")


#####################
# families
#####################

def get_families_list(datadir):
  return os.listdir(get_families_dir(datadir))

def get_family_path(datadir, family):
  return os.path.join(get_families_dir(datadir), family)

#####################
# per-family directories
#####################

def get_gene_tree_dir(datadir, family):
  return os.path.join(get_family_path(datadir, family), "gene_trees")

def get_family_misc_dir(datadir, family):
  return os.path.join(get_family_path(datadir, family), "misc")

def get_mappings_dir(datadir, family):
  return os.path.join(get_family_path(datadir, family), "mappings")

#####################
# per-family files
#####################

def build_gene_tree_path_from_run(datadir, family, run):
  return os.path.join(get_gene_tree_dir(datadir, family), run + ".geneTree.newick")

def build_gene_tree_path(datadir, subst_model, family, method):
  if (method == "true"):
    return os.path.join(get_gene_tree_dir(datadir, family), method + "." + method + ".geneTree.newick")
  else:
    return os.path.join(get_gene_tree_dir(datadir, family), method + "." + subst_model + ".geneTree.newick")

def build_misc_file_path(datadir, subst_model, family, method):
  return os.path.join(get_family_misc_dir(datadir, family), method + "." + subst_model + ".txt")

def get_results(datadir, family):
  return os.path.join(get_family_path(datadir, family), "results")

def get_alignment(datadir, family):
  return os.path.join(get_family_path(datadir, family), "alignment.msa")

def get_true_tree(datadir, family):
  return build_gene_tree_path(datadir, "true", family, "true")

def get_raxml_tree(datadir, subst_model, family):
  return build_gene_tree_path(datadir, subst_model, family, "raxml-ng")

def get_raxml_light_tree(datadir, subst_model, family):
  return build_gene_tree_path(datadir, subst_model, family, "raxml-light")

def get_raxml_multiple_trees(datadir, subst_model, family):
  return build_gene_tree_path(datadir, subst_model, family, "raxmlMultiple")

def get_phyldog_tree(datadir, subst_model, family):
  return build_gene_tree_path(datadir, subst_model, family, "phyldog")

def get_treerecs_tree(datadir, subst_model, family):
  return build_gene_tree_path(datadir, subst_model, family, "treerecs")

def get_treefix_tree(datadir, subst_model, family):
  return build_gene_tree_path(datadir, subst_model, family, "treefix")

def get_eccetera_tree(datadir, subst_model, family):
  return build_gene_tree_path(datadir, subst_model, family, "eccetera")

def get_deleterious_tree(datadir, subst_model, family):
  return build_gene_tree_path(datadir, subst_model, family, "deleterious")

def get_notung_tree(datadir, subst_model, family, threshold):
  return build_gene_tree_path(datadir, subst_model, family,  "notung" + str(threshold))

def get_ale_tree(datadir, subst_model, family, method):
  return build_gene_tree_path(datadir, subst_model, family, method)

def get_mappings(datadir, family):
  return os.path.join(get_mappings_dir(datadir, family), "mapping.link")

def get_treerecs_mappings(datadir, family):
  return os.path.join(get_mappings_dir(datadir, family), "treerecs_mapping.link")

def get_alignment_file(datadir):
  return os.path.join(datadir, "alignment.msa")

def get_raxml_best_model(datadir, subst_model, family):
  return build_misc_file_path(datadir, subst_model, family, "raxmlBestModel")

def get_raw_rf_cells_file(datadir, rooted = False):
  if (rooted):
    return os.path.join(datadir, "misc", "rooted_rf_cells.pickle")
  else:
    return os.path.join(datadir, "misc", "rf_cells.pickle")

def get_raw_kf_cells_file(datadir, rooted = False):
  return os.path.join(datadir, "misc", "kf_cells.pickle")


def get_gene_tree(datadir, subst_model, family, tree):
  gene_trees_dir = get_gene_tree_dir(datadir, family)
  lower_tree = tree.lower()
  if (lower_tree == "raxml-ng"):
    return get_raxml_tree(datadir, subst_model, family)
  elif (lower_tree == "raxmls"):
    return get_raxml_multiple_trees(datadir, subst_model, family)
  elif (lower_tree == "true"):
    return get_true_tree(datadir, family)
  elif (lower_tree == "treerecs"):
    return os.path.join(gene_trees_dir, "treerecsGeneTree.newick")
  elif (lower_tree == "phyldog"):
    return os.path.join(gene_trees_dir, "phyldogGeneTree.newick")
  elif ("notung" in lower_tree):
    return os.path.join(gene_trees_dir, lower_tree + "GeneTree.newick")
  elif ("ale" in lower_tree):
    return os.path.join(gene_trees_dir, tree + "GeneTree.newick")
  elif (lower_tree == "random"):
    res = os.path.join(gene_trees_dir, "randomGeneTree.newick")
    if (os.path.isfile(res)):
        return res
    else:
      return "__random__";
  else:
    return build_gene_tree_path(datadir, subst_model, family, lower_tree)

def get_method_name_aux(tree_path):
  tree = os.path.basename(tree_path)
  if (tree == "raxmlGeneTree.newick"):
    return "raxml-ng"
  elif (tree == "raxmlGeneTrees.newick"):
    return "raxmls"
  elif (tree == "trueGeneTree.newick"):
    return "true"
  elif (tree == "treerecsGeneTree.newick"):
    return "treerecs"
  elif (tree == "phyldogGeneTree.newick"):
    return "phyldog"
  elif (tree == "randomGeneTree.newick"):
    return "random"
  elif (tree == "randomGeneTree.newick"):
    return "random"
  else:
    return tree.replace("GeneTree.newick", "").replace(".newick", "")

def get_method_from_run(run):
  return run.split(".")[0]

def get_run_name(method, model):
  return method + "." + model

def get_run_name_from_tree(tree_file):
  split = os.path.basename(tree_file).split(".")
  return split[0] + "." + split[1]

def get_successful_runs(datadir):
  runs = []
  families_dir = os.path.join(datadir, "families")
  one_family_dir = os.path.join(families_dir, os.listdir(families_dir)[0])
  directories_to_check = []
  directories_to_check.append(os.path.join(one_family_dir, "gene_trees"))
  for directory in directories_to_check:
    for tree in os.listdir(directory):
      run_name = get_run_name_from_tree(tree)
      if (not "raxmlMultiple" in run_name):
        runs.append(run_name)
  return runs

def get_method_name(family_dir, tree_path):
  if (not tree_path.endswith(".newick")):
    return None
  method_name = get_method_name_aux(tree_path)
  if (not os.path.isfile(get_gene_tree(family_dir, method_name))):
    return None
  return method_name


def get_ran_methods(datadir):
  methods = []
  families_dir = os.path.join(datadir, "families")
  one_family_dir = os.path.join(families_dir, os.listdir(families_dir)[0])
  print(one_family_dir)
  directories_to_check = []
  directories_to_check.append(os.path.join(one_family_dir, "gene_trees"))

  for directory in directories_to_check:
    for tree in os.listdir(directory):
      method = get_method_name(one_family_dir, os.path.join(directory, tree))
      if (method != None and method != "raxmls"):
        methods.append(method)
  print(methods)
  return methods

def get_param_from_dataset_name(parameter, dataset):
  if (parameter == "species"):
    return dataset.split("_")[1][1:]
  if (parameter == "families"):
    return dataset.split("_")[2][1:]
  elif (parameter == "sites"):
    return dataset.split("_")[3][5:]
  elif (parameter == "model"):
    return dataset.split("_")[4]
  elif (parameter == "bl"):
    return dataset.split("_")[5][2:]
  elif (parameter == "dup_rate"):
    return dataset.split("_")[6][1:]
  elif (parameter == "loss_rate"):
    return dataset.split("_")[7][1:]
  elif (parameter == "transfer_rate"):
    return dataset.split("_")[8][1:]
  elif (parameter == "perturbation"):
    return dataset.split("_")[9][1:]
  elif (parameter == "population"):
    return dataset.split("_")[10][3:]
  elif (parameter == "seed"):
    return dataset.split("_")[11][4:]
  elif (parameter == "tl_ratio"):
    t = get_param_from_dataset_name("transfer_rate", dataset)
    l = get_param_from_dataset_name("loss_rate", dataset)
    return  str(float(t)/float(l))
  elif (parameter == "dl_ratio"):
    d = get_param_from_dataset_name("dup_rate", dataset)
    l = get_param_from_dataset_name("loss_rate", dataset)
    return str(float(d)/float(l))
  elif (parameter == "dt_ratio"):
    d = get_param_from_dataset_name("dup_rate", dataset)
    t = get_param_from_dataset_name("transfer_rate", dataset)
    return str(float(d)/float(t))
  elif (parameter == "av_rate"):
    d = float(get_param_from_dataset_name("dup_rate", dataset))
    l = float(get_param_from_dataset_name("loss_rate", dataset))
    t = float(get_param_from_dataset_name("transfer_rate", dataset))
    return str((d + t + l) / 2.0)
  else:
    return "invalid"

######################
# convert functions
######################
def convert_to_phyldog_species_tree(speciesTree, phyldog_species_tree):
  command = "sed s/)[nHR][0123456789a-zA-Z]*/)/g " + speciesTree 
  with open(phyldog_species_tree, "w") as output:
    subprocess.check_call(command.split(" "), stdout=output)
  print(open(speciesTree).read())
  subprocess.check_call(command.split(" "))
  print("")
  print("")

def convert_phyldog_to_treerecs_mapping(phyldog_mappings, treerecs_mappings):
  lines = open(phyldog_mappings).readlines()
  with open(treerecs_mappings, "w") as writer:
    for line in lines:
      split = line.split(":")
      species = split[0]
      genes = split[1].split(";")
      for gene in genes:
        writer.write(gene.replace("\n", "") + " " + species + "\n")

def write_phyldog_mapping(species_to_genes_dict, output_file):
  with open(output_file, "w") as writer:
    for species in species_to_genes_dict:
      writer.write(species + ":" + ";".join(species_to_genes_dict[species]) + "\n")

######################
#  Misc
######################

def get_lb_from_run(run_dir):
  lines = open(os.path.join(run_dir, "statistics.svg")).readlines()
  for line in lines:
    if (", lb = " in line):
      return float(line.split("lb = ")[1].split("<")[0])
  assert(False)
  return None


#######################
#  Directory helpers
######################

def init_top_directories(datadir):
  mkdir(datadir)
  mkdir(get_species_dir(datadir))
  mkdir(get_adjacencies_dir(datadir))
  mkdir(get_families_dir(datadir))
  mkdir(get_alignments_dir(datadir))
  mkdir(get_misc_dir(datadir))
  mkdir(get_run_dir(datadir))

def init_family_directories(datadir, family):
  mkdir(get_gene_tree_dir(datadir, family))
  mkdir(get_family_misc_dir(datadir, family))
  mkdir(get_mappings_dir(datadir, family))

def init_families_directories(datadir, families):
  for family in families:
    init_family_directories(datadir, family)


def postprocess_datadir(datadir):
  # phyldog species trees
  convert_to_phyldog_species_tree(get_species_tree(datadir), get_phyldog_species_tree(datadir)) 
  # alignments
  for family in get_families_list(datadir):
    family_dir = get_family_path(datadir, family)
    exp.relative_symlink(get_alignment(datadir, family), os.path.join(datadir, "alignments", family + ".fasta"))
    convert_phyldog_to_treerecs_mapping(get_mappings(datadir, family), get_treerecs_mappings(datadir, family)) 

def perturbate_species_tree(datadir, perturbation):
  shutil.copyfile(get_species_tree(datadir), get_true_species_tree(datadir))
  perturbate.perturbate(get_true_species_tree(datadir), get_species_tree(datadir), perturbation)




