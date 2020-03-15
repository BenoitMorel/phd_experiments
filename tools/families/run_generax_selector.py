import sys
import os
sys.path.insert(0, 'scripts/generax')
sys.path.insert(0, 'tools/family')

import launch_generax
import fam
import shutil

def select(dataset_dir, subst_model, method_names_file, cores):
  method_names = open(method_names_file).read().replace("\n", "").split(" ")
  likelihoods = []
  run_dir = "plop"
  try:
    shutil.rmtree(run_dir)
    os.makedirs(run_dir)
  except:
    pass
  for method in method_names:
    print("Evaluating  joint likelihood for " + method + " species tree")
    resultsdir = os.path.join(run_dir, method)
    os.makedirs(resultsdir)
    dataset = os.path.basename(os.path.normpath(dataset_dir))
    starting_tree = "raxml-ng"
    optimize_species = False
    radius = 5
    additional_arguments = ["--max-spr-radius", str(radius)] 
    output_dir = launch_generax.run(dataset, subst_model, "SPR", method, starting_tree, cores, additional_arguments, resultsdir, False, False)
    stats_file = os.path.join(output_dir, "generax", "stats.txt")
    ll = float(open(stats_file).readline().split()[1])
    print(ll)
    likelihoods.append((ll, method))
  likelihoods = sorted(likelihoods, reverse = True)
  best_method_with_true = likelihoods[0][1]
  best_method = best_method_with_true
  if (best_method == "true"):
    best_method = likelihoods[1][1]
  select_dest = fam.get_species_tree(dataset_dir, subst_model, "generax_select")
  select_true_dest = fam.get_species_tree(dataset_dir, subst_model, "generax_select_true")
  best_tree_with_true = fam.get_species_tree(dataset_dir, subst_model, best_method_with_true) 
  best_tree = fam.get_species_tree(dataset_dir, subst_model, best_method) 
  print(likelihoods)
  shutil.copy(best_tree_with_true, select_true_dest)
  shutil.copy(best_tree, select_dest)

if (__name__== "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax: python " + os.path.basename(__file__) + " dataset_dir subst_model method_names_file cores")
    sys.exit(1)
  dataset_dir = sys.argv[1]
  subst_model = sys.argv[2]
  method_names_file = sys.argv[3]
  cores = int(sys.argv[4])
  select(dataset_dir, subst_model, method_names_file, cores)
