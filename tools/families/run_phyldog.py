import os
import sys
import subprocess
import shutil
import time
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "trees"))
import saved_metrics
import fam
import experiments as exp
import nhx_to_newick
import sequence_model

def get_phyldog_run_name(opt_species_tree, subst_model):
  res = "Phyldog"
  if (opt_species_tree):
    res += "Species"
  return fam.get_run_name(res, subst_model)

def get_phyldog_run_dir(opt_species_tree, datadir, subst_model):
  if (opt_species_tree):
    return fam.get_run_dir(datadir, subst_model, "phyldog_species_run")
  else:
    return fam.get_run_dir(datadir, subst_model, "phyldog_run")

def clean_phyldog(opt_species_tree, datadir, subst_model):
  phyldog_run_dir = get_phyldog_run_dir(opt_species_tree, datadir, subst_model)
  dataset_name = os.path.basename(datadir)
  print("CLEANING " + dataset_name)
  for f in os.listdir(phyldog_run_dir):
    if (f.startswith("tmpPLL") and  (dataset_name.replace(".", "_") in f.replace(".", "_"))):
      os.remove(os.path.join(phyldog_run_dir, f))

def run_phyldog_on_families(datadir, subst_model, cores, opt_species_tree = False):
  output_dir = get_phyldog_run_dir(opt_species_tree, datadir, subst_model)
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  generate_options(datadir, subst_model, opt_species_tree)
  start = time.time()
  run_phyldog(datadir, subst_model, opt_species_tree, cores)
  saved_metrics.save_metrics(datadir, get_phyldog_run_name(opt_species_tree, subst_model), (time.time() - start), "runtimes") 
  saved_metrics.save_metrics(datadir, get_phyldog_run_name(opt_species_tree, subst_model), (time.time() - start), "seqtimes") 
  extract_phyldog(datadir, subst_model, opt_species_tree)
  clean_phyldog(opt_species_tree, datadir, subst_model)


def add_starting_tree(option_file, tree_path):
  lines = open(option_file).readlines()
  with open(option_file, "w") as writer:
    for line in lines:
      if (line.startswith("init.gene.tree")):
        continue
      if (line.startswith("gene.tree.file")):
          continue
      writer.write(line)
    writer.write("\n")
    writer.write("init.gene.tree=user\n")
    writer.write("gene.tree.file=" + tree_path + "\n")
    writer.write("use.quality.filters=0\n")

def generate_options(datadir, subst_model, opt_species_tree):
  phyldog_run_dir = get_phyldog_run_dir(opt_species_tree, datadir, subst_model)
  prepare_input = os.path.join(phyldog_run_dir, "prepare_input.txt")
  datadir = os.path.abspath(datadir)
  phyldog_run_dir = os.path.abspath(phyldog_run_dir)
  mappings_dir = os.path.join(phyldog_run_dir, "mappings")
  exp.reset_dir(mappings_dir)
  all_raxml_trees_dir = os.path.join(phyldog_run_dir, "all_raxml_trees")
  exp.reset_dir(all_raxml_trees_dir)
  for family in fam.get_families_list(datadir):
    phyldog_mapping = fam.get_mappings(datadir, family)
    new_phyldog_mapping = os.path.join(mappings_dir, family + ".link")
    exp.relative_symlink(phyldog_mapping, new_phyldog_mapping)
    old_raxml_tree = fam.get_raxml_tree(datadir, subst_model, family)
    new_raxml_tree = os.path.join(all_raxml_trees_dir, family + ".newick")
    exp.relative_symlink(old_raxml_tree, new_raxml_tree) 
  with open(prepare_input, "w") as writer:
    writer.write(os.path.join(datadir, "alignments") + "\n")
    if (sequence_model.is_dna(subst_model)):
      writer.write("1\n" )
    else:
      writer.write("0\n")
    model = sequence_model.get_phyldog_model(subst_model)
    if (model == "LG"):
      model = "LG08"
    if (model == "WAG"):
      model = "WAG01"
    writer.write(model + "\n")
    writer.write(str(sequence_model.get_gamma_rates(subst_model)) + "\n")
    writer.write("FASTA\n")
    writer.write(os.path.join(phyldog_run_dir, "mappings") + "\n")
    writer.write(os.path.join(phyldog_run_dir, "options") + "\n")
    writer.write(os.path.join(phyldog_run_dir, "results") + "\n")
    writer.write("yes" + "\n")
    writer.write(fam.get_phyldog_species_tree(datadir) + "\n")
    if (opt_species_tree):
      writer.write("yes" + "\n") #opt species tree
    else:
      writer.write("no" + "\n") #opt species tree
    writer.write("yes" + "\n") # opt dup loss
    writer.write("average" + "\n") # branchwise DL opt
    writer.write("no" + "\n") # same number of genes ?
    writer.write("yes" + "\n") #opt gene trees
    writer.write("48" + "\n") # max time (hours)
  prepare_data_script = os.path.join(exp.tools_root, "families", "prepareData.py")
  logs = open(os.path.join(phyldog_run_dir, "options_logs.txt"), "w")
  
  subprocess.check_call(["/bin/bash", "-c", "python " + prepare_data_script + " < " + prepare_input], stdout = logs)
  time.sleep(1)
  for family in fam.get_families_list(datadir):
    option_file = os.path.join(phyldog_run_dir, "options", family + ".opt")
    add_starting_tree(option_file, fam.get_raxml_tree(datadir, subst_model, family))


   
def run_phyldog(datadir, subst_model, opt_species_tree, cores):
  cwd = os.getcwd()
  try:
    families_number = len(os.listdir(os.path.join(datadir, "families")))
    cores = str(min(families_number, int(cores)))
    phyldog_run_dir = get_phyldog_run_dir(opt_species_tree, datadir, subst_model)
    command = []
    command.append("mpirun")
    command.append("-n")
    command.append(str(cores))
    command.append(exp.phyldog_exec)
    command.append("param=" + os.path.abspath(os.path.join(phyldog_run_dir, "options", "GeneralOptions.txt")))
    print("EXECUTE PHYLDOG")
    print(" ".join(command))
    logs = open(os.path.join(phyldog_run_dir, "logs.txt"), "w")
    os.chdir(phyldog_run_dir)
    start = time.time()
    subprocess.check_call(command, stdout = logs, stderr = logs)
    saved_metrics.save_metrics(datadir, "Phyldog", (time.time() - start), "runtimes") 
    print("end EXECUTE PHYLDOG")
  finally:
    os.chdir(cwd)


def extract_phyldog(datadir, subst_model, opt_species_tree):
  results_dir = os.path.join(get_phyldog_run_dir(opt_species_tree, datadir, subst_model), "results")
  families_dir = os.path.join(datadir, "families")
  for family in os.listdir(families_dir):
    phyldog_tree = os.path.join(results_dir, family + ".ReconciledTree")
    try:
      nhx_to_newick.nhx_to_newick(phyldog_tree, fam.get_phyldog_tree(datadir, subst_model, family))
    except:
      print("Phyldog failed to infer tree " + phyldog_tree)
      shutil.copy(fam.get_raxml_tree(datadir, family), fam.get_phyldog_tree(datadir, subst_model, family))
  if (opt_species_tree):
    species_tree = os.path.join(results_dir, "OutputSpeciesTree_ConsensusDuplications.tree")
    new_species_tree = fam.get_species_tree(datadir, None, get_phyldog_run_name(opt_species_tree, subst_model).lower())
    print("Saving phyldog tree in " + new_species_tree)
    shutil.copy(species_tree, new_species_tree)

if (__name__== "__main__"):
  max_args_number = 5
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_phyldog.py datadir subst_model optimize_species_tree cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  opt_species_tree = int(sys.argv[3]) != 0
  cores = int(sys.argv[4])
  run_phyldog_on_families(datadir, subst_model, cores, opt_species_tree)

