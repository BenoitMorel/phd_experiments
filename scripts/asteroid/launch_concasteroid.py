mport sys
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
import run_fastme


def get_distance_matrix_file(datadir, subst_model, output_dir):
  f = os.path.join(output_dir, "distance_matrices.txt")
  with open(f, "w") as writer:
    for family in fam.get_families_list(datadir):
      fastme_matrix = fam.get_fastme_distances(datadir, family, subst_model)
      if (os.path.isfile(fastme_matrix)):
        writer.write(fastme_matrix)
        writer.write("\n")
  return f


def get_mapping_file(datadir, output_dir, additional_arguments):
  if ("--single" in additional_arguments):
    additional_arguments.remove("--single")
    return None
  mapping_file = os.path.join(output_dir, "mapping_file.txt")
  m = get_dico.get_species_to_genes(datadir)
  #print(m)
  print(mapping_file)
  get_dico.export_species_to_genes_to_file(m, mapping_file)
  return mapping_file


def get_prefix(output_dir):
  return  os.path.join(output_dir, "asteroid")

def get_inferred_tree(output_dir):
  return  get_prefix(output_dir) + ".bestTree.newick"

def get_asteroid_command(distance_matrices, mapping_file, additional_arguments, output_dir, cores):
    executable = exp.concasteroid_exec
    prefix = get_prefix(output_dir)
    command = []
    command.append("mpirun")
    command.append("-np")
    command.append(str(cores))
    command.append(executable)
    command.append("-p")
    command.append(prefix)
    command.append("-i")
    command.append(distance_matrices)
    if (mapping_file != None):
      command.append("-m")
      command.append(mapping_file)
    command.extend(additional_arguments)
    return " ".join(command)


def run_concasteroid(datadir, cores, distance_matrices, mapping_file, additional_arguments, output_dir):
  command = get_asteroid_command(distance_matrices, mapping_file, additional_arguments, output_dir, cores)
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

def run(datadir, subst_model, is_dna, cores, additional_arguments):
  run_name = exp.getAndDelete("--run", additional_arguments, None) 
  if (run_name == None):
    noCorrection = "-n" in additional_arguments 
    noCorrection = noCorrection or "--no-correction" in additional_arguments
    run_name = "concasteroid-"
    if (noCorrection):
      run_name = "concfastme-"
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
    run_name += "." + subst_model
  
  output_dir = fam.get_run_dir(datadir, subst_model, "concasteroid_run")
  try:
    os.makedirs(output_dir)
  except:
    print("Failed to create " + output_dir)
    pass
  if (not exp.checkAndDelete("--skip", additional_arguments)):
    run_fastme.run_fastme_on_families(datadir, subst_model, is_dna, cores)
  start = time.time()
  mapping_file = get_mapping_file(datadir, output_dir, additional_arguments)
  distance_matrices = get_distance_matrix_file(datadir, subst_model, output_dir)
  run_concasteroid(datadir, cores, distance_matrices, mapping_file, additional_arguments, output_dir)
  saved_metrics.save_metrics(datadir, run_name, (time.time() - start), "runtimes") 
  extract_trees(datadir, output_dir, run_name, subst_model)
  analyze_species_results(datadir, output_dir)
  print("Output in " + output_dir)

  

if (__name__ == "__main__"): 
    
  min_args_number = 5
  if (len(sys.argv) < min_args_number):
    print("Syntax error: python " + os.path.basename(__file__) + "  datadir subst_model subst_model is_dna cores [additional paremeters].\n")
    sys.exit(1)
  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  is_dna = sys.argv[3] != "0"
  cores = int(sys.argv[4])
  dataset = os.path.basename(os.path.normpath(datadir))
  if (dataset == datadir):
    datadir = fam.get_datadir(dataset)
  additional_arguments = sys.argv[min_args_number:]
  run(datadir, subst_model, is_dna, cores, additional_arguments)



