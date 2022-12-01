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




def av_rf(rf_cell):
  return float(rf_cell[0]) / float(rf_cell[1])

def analyze_species_results(datadir, resultsdir):
  true_species_tree = Tree(fam.get_species_tree(datadir), format = 1)
  starting_species_tree = Tree(os.path.join(resultsdir, "alegenerax", "species_trees", "starting_species_tree.newick"), format = 1)
  inferred_species_tree = Tree(os.path.join(resultsdir, "alegenerax", "species_trees", "inferred_species_tree.newick"), format = 1)
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

def get_starting_gene_tree_path(datadir, subst_model, family, starting_gene_tree):
  return fam.build_gene_tree_path(datadir, subst_model, family, starting_gene_tree)

def do_not_opt_rates(additional_arguments):
  arg = "--dtl-rates-opt"
  if (not arg in additional_arguments):
    return False
  pos = additional_arguments.index(arg)
  return (additional_arguments[pos + 1] == "NONE")

def export_gene_trees(datadir, output_dir, run_name):
  reconciliations = os.path.join(output_dir, "alegenerax", "reconciliations")
  for family in fam.get_families_list(datadir):
    source = os.path.join(reconciliations, family + ".newick")
    dest = fam.build_gene_tree_path_from_run(datadir, family, run_name)
    shutil.copy(source, dest)


def build_alegenerax_families_file(datadir, starting_gene_tree, subst_model, output):
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

def get_alegenerax_command(alegenerax_families_file, starting_species_tree, additional_arguments, output_dir, mode, cores):
    executable = exp.genetegrator_exec
    alegenerax_output = os.path.join(output_dir, "alegenerax")
    command = []
    command.append("mpirun")
    command.append("-np")
    command.append(str(cores))
    command.append(executable)
    command.append("-f")
    command.append(alegenerax_families_file)
    command.append("-s")
    command.append(starting_species_tree)
    command.append("-p")
    command.append(alegenerax_output)
    command.append("--species-search")
    command.append("EVAL")
    command.append("--gene-tree-samples")
    command.append("1")
    command.extend(additional_arguments)
    return " ".join(command)

def run_alegenerax(datadir, starting_species_tree, alegenerax_families_file, mode, cores, additional_arguments, resultsdir):
  species_tree = fam.get_species_tree(datadir)
  command = get_alegenerax_command(alegenerax_families_file, starting_species_tree, additional_arguments, resultsdir, mode, cores)
  print("Running:")
  print(command)
  subprocess.check_call(command.split(" "), stdout = sys.stdout)

def extract_species_tree(datadir, results_family_dir, run_name, subst_model):
  src = os.path.join(results_family_dir, "species_trees", "inferred_species_tree.newick")
  dest = fam.get_species_tree(datadir, None, run_name)
  print("extracted tree " + dest)
  shutil.copyfile(src, dest)



def run(datadir, subst_model, starting_species_tree, starting_gene_tree, cores, additional_arguments, resultsdir, do_analyze = True, do_extract = True):
  run_name = exp.getAndDelete("--run", additional_arguments, None) 
  
  if (run_name == None):
    run_name = "alegenerax" 
    tc = exp.getArg("--transfer-constraint", additional_arguments, "NONE")
    if (tc == "NONE"):
      run_name += "-tcnone"
    if (tc == "SOFTDATED"):
      run_name += "-tcsoft"
    rec_model = exp.getArg("--rec-model", additional_arguments, "UndatedDTL")
    if (rec_model != "UndatedDTL"):
      run_name += rec_model
    gamma = exp.getArg("--gamma-categories", additional_arguments, "1")
    if (gamma != "1"):
      run_name += "_cat" + gamma
    run_name += "_" + starting_gene_tree
    run_name += "." + subst_model
    
  arg_analyze = exp.getAndDelete("--analyze", additional_arguments, "yes")
  do_analyze = do_analyze and (arg_analyze == "yes")
  print("Run name " + run_name)
  sys.stdout.flush()
  alegenerax_families_file = os.path.join(resultsdir, "families.txt")
  build_alegenerax_families_file(datadir, starting_gene_tree, subst_model, alegenerax_families_file)
  start = time.time()
  species_tree = fam.get_species_tree(datadir, subst_model, starting_species_tree) 
  mode = ""
  run_alegenerax(datadir, species_tree, alegenerax_families_file, mode, cores, additional_arguments, resultsdir)
  saved_metrics.save_metrics(datadir, run_name, (time.time() - start), "runtimes") 
  saved_metrics.save_metrics(datadir, run_name, (time.time() - start), "seqtimes") 
  if (do_extract):
    export_gene_trees(datadir, resultsdir, run_name)
    extract_species_tree(datadir, os.path.join(resultsdir, "alegenerax"), run_name, subst_model)
  analyze_species_results(datadir, resultsdir)
  if (do_analyze):
    print("Analyzing gene trees...")
    fast_rf_cells.analyze(datadir, "all", cores, run_name)
  print("Output in " + resultsdir)

def launch(datadir, subst_model, starting_species_tree, starting_gene_tree, cluster, cores, additional_arguments):
  command = [exp.python()]
  command.extend(sys.argv)
  command.append("--exprun")
  dataset = os.path.basename(datadir)
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
    print("Syntax error: python " + os.path.basename(__file__) + "  datadir species_tree gene_trees subst_model cluster cores [additional paremeters]")
    sys.exit(1)

  datadir = os.path.normpath(sys.argv[1])
  starting_species_tree = sys.argv[2]
  starting_gene_tree = sys.argv[3]
  subst_model = sys.argv[4]
  cluster = sys.argv[5]
  cores = int(sys.argv[6])
  additional_arguments = sys.argv[min_args_number:]


  if (is_run):
    run(datadir, subst_model, starting_species_tree, starting_gene_tree, cores, additional_arguments, resultsdir)
  else:
    launch(datadir, subst_model, starting_species_tree, starting_gene_tree, cluster, cores, additional_arguments)




