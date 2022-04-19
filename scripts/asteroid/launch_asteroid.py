import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/mappings')
import saved_metrics
import rf_cells
import experiments as exp
import shutil
import time
import fam
import sequence_model
import fast_rf_cells
from ete3 import Tree
import get_dico





def get_mapping_file(datadir, additional_arguments):
  if ("--single" in additional_arguments):
    additional_arguments.remove("--single")
    return None
  mapping_file = os.path.join(output_dir, "mapping_file.txt")
  m = get_dico.get_species_to_genes(datadir)
  #print(m)
  print(mapping_file)
  get_dico.export_species_to_genes_to_file(m, mapping_file)
  return mapping_file

def get_gene_tree_file(datadir,gene_tree_method, subst_model):
  gene_tree_file = os.path.join(output_dir, "gene_trees.txt")
  failures = 0
  with open(gene_tree_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      try:
        gene_tree = fam.build_gene_tree_path(datadir, subst_model, family, gene_tree_method)
        writer.write(open(gene_tree).read())
        writer.write("\n")
      except:
        #print("Can't load gene tree for family " + family)
        failures = failures + 1
  if (failures > 9):
    print("Failed to load " + str(failures) + " gene trees")
  return gene_tree_file

def get_prefix(output_dir):
  return  os.path.join(output_dir, "asteroid")

def get_inferred_tree(output_dir):
  return  get_prefix(output_dir) + ".bestTree.newick"

def get_asteroid_command(gene_tree_file, mapping_file, additional_arguments, output_dir, cores):
    executable = exp.asteroid_exec
    prefix = get_prefix(output_dir)
    command = []
    command.append("mpirun")
    command.append("-np")
    command.append(str(cores))
    command.append(executable)
    command.append("-p")
    command.append(prefix)
    command.append("-i")
    command.append(gene_tree_file)
    if (mapping_file != None):
      command.append("-m")
      command.append(mapping_file)
    command.extend(additional_arguments)
    return " ".join(command)


def run_asteroid(datadir, cores, gene_tree_file, mapping_file, additional_arguments, output_dir):
  command = get_asteroid_command(gene_tree_file, mapping_file, additional_arguments, output_dir, cores)
  print("Running:")
  print(command)
  try:
    subprocess.check_call(command.split(" "), stdout = sys.stdout)
  except:
    print("Command: ")
    print(command)
    print("failed!")
    sys.exit(0)
  #print("I abort, remove this line!")
  #sys.exit(0)

def extract_trees(datadir, output_dir, run_name, subst_model):
  src = get_inferred_tree(output_dir)
  print(src)
  dest = fam.get_species_tree(datadir, None, run_name)
  shutil.copyfile(src, dest)
  print("extracted tree " + dest)

def av_rf(rf_cell):
  if (float(rf_cell[1]) == 0.0):
    return 10000000.0
  return float(rf_cell[0]) / float(rf_cell[1])


def analyze_species_results(datadir, output_dir):
  true_species_tree = Tree(fam.get_species_tree(datadir), format = 1)
  inferred_species_tree = Tree(get_inferred_tree(output_dir), format = 1)
  inferred_unrooted_rf = true_species_tree.robinson_foulds(inferred_species_tree, unrooted_trees = True, correct_by_polytomy_size = True)
  print("Inferred species tree:")
  print("  Unrooted RF: " + str(av_rf(inferred_unrooted_rf)))

def run(datadir, gene_tree_method, subst_model, cores, additional_arguments, output_dir):
  run_name = exp.getAndDelete("--run", additional_arguments, None) 
  if (run_name == None):
    noCorrection = "-n" in additional_arguments 
    noCorrection = noCorrection or "--no-correction" in additional_arguments
    run_name = "asteroid-"
    if (noCorrection):
      run_name = "asteroidastrid-"
    if ("--stepwise" in additional_arguments):
      run_name += "stepwise"
    if (int(cores) > 1):
      run_name += "cores" + str(cores) + "-"
    if ("-r" in additional_arguments):
      random_starting_trees = exp.getArg("-r", additional_arguments, "1")
      run_name += "r" + str(random_starting_trees) + str("-")
    if ("--min-bl" in additional_arguments):
      minbl = exp.getArg("--min-bl", additional_arguments, "1")
      run_name += "minbl" + str(minbl) + "-"
    if ("--seed" in additional_arguments):
      seed = exp.getArg("--seed", additional_arguments, "1")
      run_name += "seed" + str(seed) + "-"
    if ("-b" in additional_arguments):
      bs = exp.getArg("-b", additional_arguments, "1")
      run_name += "b" + str(bs) + "-"
    run_name += gene_tree_method
    run_name += "." + subst_model
   
  start = time.time()
  gene_tree_file = get_gene_tree_file(datadir,gene_tree_method, subst_model)
  mapping_file = get_mapping_file(datadir, additional_arguments)
  run_asteroid(datadir, cores, gene_tree_file, mapping_file, additional_arguments, output_dir)
  saved_metrics.save_metrics(datadir, run_name, (time.time() - start), "runtimes") 
  extract_trees(datadir, output_dir, run_name, subst_model)
  analyze_species_results(datadir, output_dir)
  print("Output in " + output_dir)

def launch(datadir, gene_tree_method, subst_model, cluster, cores, additional_arguments, output_dir):
  dataset = os.path.basename(os.path.normpath(datadir))
  command = [exp.python()]
  command.extend(sys.argv)
  command.append("--exprun")
  output_dir = os.path.join("Asteroid", dataset, gene_tree_method, "run")
  output_dir = exp.create_result_dir(output_dir, additional_arguments)
  submit_path = os.path.join(output_dir, "submit.sh")
  command.append(output_dir)
  exp.submit(submit_path, " ".join(command), cores, cluster) 
  

if (__name__ == "__main__"): 
  is_run = ("--exprun" in sys.argv)
  output_dir = ""
  if (is_run):
    output_dir = sys.argv[-1]
    sys.argv = sys.argv[:-2]
    
  min_args_number = 6
  if (len(sys.argv) < min_args_number):
    print("Syntax error: python " + os.path.basename(__file__) + "  datadir gene_tree_method subst_model cluster cores [additional paremeters].\n")
    sys.exit(1)

  datadir = sys.argv[1]
  gene_tree_method = sys.argv[2]
  subst_model = sys.argv[3]
  cluster = sys.argv[4]
  cores = int(sys.argv[5])
  dataset = os.path.basename(os.path.normpath(datadir))
  if (dataset == datadir):
    datadir = fam.get_datadir(dataset)
  additional_arguments = sys.argv[min_args_number:]

  if (is_run):
    run(datadir, gene_tree_method, subst_model, cores, additional_arguments, output_dir)
  else:
    launch(datadir, gene_tree_method, subst_model, cluster, cores, additional_arguments, output_dir)


