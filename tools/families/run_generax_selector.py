import sys
import os
import time
sys.path.insert(0, 'scripts/generax')
sys.path.insert(0, 'tools/family')

import launch_generax
import fam
import shutil
import saved_metrics

def need_reroot(method_name):
  if (method_name == "true"):
    return False
  if ("speciesrax" in method_name or "generax" in method_name):
    return False
  return True

def select(dataset_dir, subst_model, method_names_file, cores):
  method_names = open(method_names_file).read().replace("\n", "").split(" ")
  likelihoods = []
  run_dir = fam.get_run_dir(dataset_dir, subst_model, "generaxselect_run")
  try:
    shutil.rmtree(run_dir)
  except:
    pass
  os.makedirs(run_dir)
  start = time.time()
  for method in method_names:
    print("Evaluating  joint likelihood for " + method + " species tree")
    resultsdir = os.path.join(run_dir, method)
    os.makedirs(resultsdir)
    dataset = os.path.basename(os.path.normpath(dataset_dir))
    starting_tree = "raxml-ng"
    optimize_species = False
    radius = 5
    additional_arguments = ["--max-spr-radius", str(radius)] 
    if (need_reroot(method)):
      additional_arguments.append("--reroot-species-tree")
    output_dir = launch_generax.run(dataset, subst_model, "SPR", method, starting_tree, cores, additional_arguments, resultsdir, False, False)
    stats_file = os.path.join(output_dir, "generax", "stats.txt")
    ll = float(open(stats_file).readline().split()[1])
    tree = os.path.join(output_dir, "generax", "inferred_species_tree.newick")
    print(ll)
    likelihoods.append((ll, method, tree))
  time1 = (time.time() - start)
  saved_metrics.save_metrics(dataset_dir, fam.get_run_name("generax-select", subst_model), time1, "runtimes") 
  likelihoods = sorted(likelihoods, reverse = True)
  best_tuple_with_true = likelihoods[0]
  best_tuple = likelihoods[0]
  if (best_tuple[1] == "true"):
    best_tuple = likelihoods[1]
  select_dest = fam.get_species_tree(dataset_dir, subst_model, "generax-select")
  select_true_dest = fam.get_species_tree(dataset_dir, subst_model, "generax-select-true")
  for t in likelihoods:
    print(str(t[0]) + " " + str(t[1]))
  shutil.copy(best_tuple_with_true[2], select_true_dest)
  shutil.copy(best_tuple[2], select_dest)

if (__name__== "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax: python " + os.path.basename(__file__) + " dataset_dir subst_model method_names_file cores")
    sys.exit(1)
  dataset_dir = sys.argv[1]
  subst_model = sys.argv[2]
  method_names_file = sys.argv[3]
  cores = int(sys.argv[4])
  select(dataset_dir, subst_model, method_names_file, cores)
