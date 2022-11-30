import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import saved_metrics
import rf_cells
import experiments as exp
import shutil
import time
import fam
import sequence_model
import fast_rf_cells
from ete3 import Tree

def get_possible_strategies():
  return ["SPR", "EVAL", "SKIP"]





def has_multiple_sample(starting_gene_tree):
  return "ale" in starting_gene_tree.lower() or "multiple" in starting_gene_tree.lower()

def get_starting_gene_tree_path(datadir, subst_model, family, starting_gene_tree):
  #if (has_multiple_sample(starting_gene_tree)):
    #return os.path.join(fam.get_family_misc_dir(datadir, family), starting_gene_tree + "." + subst_model + "_onesample.geneTree")
  #else:
  return fam.build_gene_tree_path(datadir, subst_model, family, starting_gene_tree)

# GeneRax does not accept tree files with multiple trees
def sample_one_starting_gene_tree(datadir, subst_model, starting_gene_tree):
  for family in fam.get_families_list(datadir):
    input_tree = fam.build_gene_tree_path(datadir, subst_model, family, starting_gene_tree)
    output_tree = get_starting_gene_tree_path(datadir, subst_model, family, starting_gene_tree)
  #  tree = open(input_tree, "r").readline()
  #  open(output_tree, "w").write(tree)

def build_generax_families_file(datadir, starting_gene_tree, subst_model, output):
  if (has_multiple_sample(starting_gene_tree)):
    sample_one_starting_gene_tree(datadir, subst_model, starting_gene_tree)
  families_dir = os.path.join(datadir, "families")
  with open(output, "w") as writer:
    writer.write("[FAMILIES]\n")
    plop = 0
    print("starting gene tree " + starting_gene_tree)
    for family in os.listdir(families_dir):
      
      family_path = os.path.join(families_dir, family)
      writer.write("- " + family + "\n")
      gene_tree = get_starting_gene_tree_path(datadir, subst_model, family, starting_gene_tree)
      if (starting_gene_tree == "random"):
        gene_tree = "__random__"
      writer.write("starting_gene_tree = " + gene_tree + "\n")
      writer.write("alignment = " + fam.get_alignment_file(family_path) + "\n")
      mapping_file = fam.get_mappings(datadir, family)
      if (os.path.isfile(mapping_file)):
        writer.write("mapping = " + fam.get_mappings(datadir, family) + "\n")
      raxml_model = ""
      if (starting_gene_tree != "random" and starting_gene_tree != "true"):
        raxml_model = fam.get_raxml_best_model(datadir, subst_model, family)
      if (os.path.isfile(raxml_model)):
        writer.write("subst_model = " + raxml_model + "\n")
      else:
        writer.write("subst_model = " + sequence_model.get_raxml_model(subst_model) + "\n")

def get_generax_command(generax_families_file, starting_species_tree, additional_arguments, output_dir, mode, cores):
    executable = exp.generax_exec
    if (mode == "gprof"):
      executable = exp.generax_gprof_exec
    elif (mode == "scalasca"):
      executable = exp.generax_scalasca_exec
    generax_output = os.path.join(output_dir, "generax")
    command = []
    command.append("mpirun")
    command.append("-np")
    command.append(str(cores))
    command.append(executable)
    command.append("-f")
    command.append(generax_families_file)
    command.append("-s")
    command.append(starting_species_tree)
    command.append("--si-strategy")
    command.append("HYBRID")
    command.append("--strategy")
    command.append("SKIP")
    #command.append("--do-not-reconcile")
    #command.append("--si-estimate-bl")
    #command.append("--si-quartet-support")
    #command.append("--si-eqpic-radius")
    #command.append("3")
    command.append("-p")
    command.append(generax_output)
    command.extend(additional_arguments)
    print("COmmand " + " ".join(command))
    return " ".join(command)

def run_generax(datadir, starting_species_tree, generax_families_file, mode, cores, additional_arguments, resultsdir):
  command = get_generax_command(generax_families_file, starting_species_tree, additional_arguments, resultsdir, mode, cores)
  subprocess.check_call(command.split(" "), stdout = sys.stdout)


def get_mode_from_additional_arguments(additional_arguments):
  mode = "normal"
  if ("--scalasca" in additional_arguments):
    mode = "scalasca"
    additional_arguments.remove("--scalasca")
  elif ("--gprof" in additional_arguments):
    mode = "gprof"
    additional_arguments.remove("--gprof")
  return mode


def extract_trees(datadir, results_family_dir, run_name, subst_model):
  src = os.path.join(results_family_dir, "species_trees", "inferred_species_tree.newick")
  dest = fam.get_species_tree(datadir, None, run_name)
  print("extracted tree " + dest)
  shutil.copyfile(src, dest)
  #
  return
  results_dir = os.path.join(results_family_dir, "results")
  for family in fam.get_families_list(datadir):
    source = os.path.join(results_dir, family, "geneTree.newick")
    dest = fam.build_gene_tree_path_from_run(datadir, family, run_name)
    try:
      shutil.copy(source, dest)
    except:
      pass

def av_rf(rf_cell):
  return float(rf_cell[0]) / float(rf_cell[1])


def analyze_species_results(datadir, resultsdir):
  true_species_tree = Tree(fam.get_species_tree(datadir), format = 1)
  starting_species_tree = Tree(os.path.join(resultsdir, "generax", "species_trees", "starting_species_tree.newick"), format = 1)
  inferred_species_tree = Tree(os.path.join(resultsdir, "generax", "species_trees", "inferred_species_tree.newick"), format = 1)
  starting_rooted_rf = -1
  inferred_rooted_rf = -1
  try:
    starting_rooted_rf = true_species_tree.robinson_foulds(starting_species_tree, unrooted_trees = False, correct_by_polytomy_size = True)
    inferred_rooted_rf = true_species_tree.robinson_foulds(inferred_species_tree, unrooted_trees = False, correct_by_polytomy_size = True)
  except:
    print("can't computed rooted distance")
  starting_unrooted_rf = true_species_tree.robinson_foulds(starting_species_tree, unrooted_trees = True, correct_by_polytomy_size = True)
  inferred_unrooted_rf = true_species_tree.robinson_foulds(inferred_species_tree, unrooted_trees = True, correct_by_polytomy_size = True)
  print("Starting species tree:")
  print("  Rooted RF: " + str(av_rf(starting_rooted_rf)))
  print("  Unrooted RF: " + str(av_rf(starting_unrooted_rf)))
  print("Inferred species tree:")
  print("  Rooted RF: " + str(av_rf(inferred_rooted_rf)))
  print("  Unrooted RF: " + str(av_rf(inferred_unrooted_rf)))

def cleanup(resultsdir):
  reconciliations = os.path.join(resultsdir, "generax", "reconciliations")
  try:
    shutil.rmtree(reconciliations)
  except:
    pass

def do_not_opt_rates(additional_arguments):
  arg = "--dtl-rates-opt"
  if (not arg in additional_arguments):
    return False
  pos = additional_arguments.index(arg)
  return (additional_arguments[pos + 1] == "NONE")

def run(dataset, subst_model, starting_species_tree, starting_gene_tree, cores, additional_arguments, resultsdir, do_analyze = True, do_extract = True):
  run_name = exp.getAndDelete("--run", additional_arguments, None) 
  if (run_name == None):
    run_name = "generax-" + starting_species_tree
    rec_model = exp.getArg("--rec-model", additional_arguments, "UndatedDTL")
    if (starting_species_tree == "random"):
      run_name += exp.getArg("--seed", additional_arguments, "noseed")
    if (rec_model != "UndatedDTL"):
      run_name += rec_model
    if ("--si-spr-radius" in additional_arguments):
      radius = exp.getArg("--si-spr-radius", additional_arguments, "1")
      run_name += "-radius" + radius
    if ("--prune-species-tree" in additional_arguments):
      run_name += "-prune"
    if ("--no-dup" in additional_arguments):
      run_name += "-nodup"
    if ("--per-family-rates" in additional_arguments):
      run_name += "-fam"
    if (do_not_opt_rates(additional_arguments)):
      run_name += "-fixed"
    if ("--si-constrained-search" in additional_arguments):
      run_name += "-constr"
    if ("--unrooted-gene-tree" in additional_arguments):
      run_name += "-unrooted"
    run_name += "_" + starting_gene_tree
    run_name += "." + subst_model
   
  print(run_name)
  arg_analyze = exp.getAndDelete("--analyze", additional_arguments, "yes")
  do_analyze = do_analyze and (arg_analyze == "no")
  print("Run name " + run_name)
  sys.stdout.flush()
  mode = get_mode_from_additional_arguments(additional_arguments)
  datadir = fam.get_datadir(dataset)
  generax_families_file = os.path.join(resultsdir, "families.txt")
  build_generax_families_file(datadir, starting_gene_tree, subst_model, generax_families_file)
  start = time.time()
  species_tree = fam.get_species_tree(datadir, subst_model, starting_species_tree)  
  if (starting_species_tree == "MiniNJ"):
    species_tree = "MiniNJ"
  print("coucou")
  print(species_tree)
  run_generax(datadir, species_tree, generax_families_file, mode, cores, additional_arguments, resultsdir)
  saved_metrics.save_metrics(datadir, run_name, (time.time() - start), "runtimes") 
  saved_metrics.save_metrics(datadir, run_name, (time.time() - start), "seqtimes") 
  if (do_extract):
    print("DO EXTRACT")
    extract_trees(datadir, os.path.join(resultsdir, "generax"), run_name, subst_model)
  analyze_species_results(datadir, resultsdir)
  try:
    if (do_analyze):
      fast_rf_cells.analyze(datadir, "all", cores, run_name)
  except:
    print("Analyze failed!!!!")
  #cleanup(resultsdir)
  print("Output in " + resultsdir)

def launch(dataset, subst_model, starting_species_tree, starting_gene_tree, cluster, cores, additional_arguments):
  command = [exp.python()]
  command.extend(sys.argv)
  command.append("--exprun")
  resultsdir = os.path.join("GeneRax", dataset, starting_species_tree + "_start_" + starting_gene_tree, "run")
  resultsdir = exp.create_result_dir(resultsdir, additional_arguments)
  submit_path = os.path.join(resultsdir, "submit.sh")
  command.append(resultsdir)
  exp.submit(submit_path, " ".join(command), cores, cluster) 
  

if (__name__ == "__main__"): 
  is_run = ("--exprun" in sys.argv)
  resultsdir = ""
  if (is_run):
    resultsdir = sys.argv[-1]
    sys.argv = sys.argv[:-2]
    
  min_args_number = 7
  if (len(sys.argv) < min_args_number):
    print("Syntax error: python " + os.path.basename(__file__) + "  dataset gene_tree subst_model starting_species_tree cluster cores [additional paremeters]. ")
    sys.exit(1)

  dataset = sys.argv[1]
  starting_gene_tree = sys.argv[2]
  subst_model = sys.argv[3]
  starting_species_tree = sys.argv[4]
  cluster = sys.argv[5]
  cores = int(sys.argv[6])
  dataset = os.path.basename(os.path.normpath(dataset))
  additional_arguments = sys.argv[min_args_number:]

  if (starting_gene_tree == "raxml"):
    print("use raxml-ng instead of raxml please")
    exit(1)

  if (is_run):
    run(dataset, subst_model, starting_species_tree, starting_gene_tree, cores, additional_arguments, resultsdir)
  else:
    launch(dataset, subst_model, starting_species_tree, starting_gene_tree, cluster, cores, additional_arguments)



