import os
import sys
import subprocess
import shutil
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/raxml/')
import experiments as exp
import raxml_get_tca_score as tca
import time
import saved_metrics
import run_raxml_supportvalues as raxml
import sequence_model

def run_pargenes(datadir, pargenes_dir, subst_model, rand_trees, parsi_trees, bs_trees, cores, arguments):
  raxml_command = ""
  run_modeltest = (subst_model == "bestAA" or subst_model == "bestNT")
  if (not run_modeltest):
    raxml_command +="--model " + sequence_model.get_raxml_model(subst_model) + " --blopt nr_safe"
  command = []
  command.append(exp.python())
  if (cores > 4):
    command.append(exp.pargenes_script)
  else:
    command.append(exp.pargenes_script_laptop)
  command.append("-a")
  command.append(os.path.join(datadir, "alignments"))
  command.append("-b")
  command.append(str(bs_trees))
  command.append("-o")
  command.append(pargenes_dir)
  command.append("-c")
  command.append(str(cores))
  command.append("-s")
  command.append(str(rand_trees))
  if (parsi_trees > 0):
    command.append("-p")
    command.append(str(parsi_trees))
  if (len(raxml_command) > 0):
    command.append("-R")
    command.append(raxml_command)
  if (run_modeltest):
    command.append("-m")
    if (subst_model == "bestAA"):
      command.append("-d")
      command.append("aa")
  command.append("--continue")
  command.extend(arguments)
  try:
    subprocess.check_call(command, stdout = sys.stdout)
  except:
    command[0] = exp.python()
    print(" ".join(command))
    subprocess.check_call(command, stdout = sys.stdout)


def gather_likelihoods(pargenes_dir, output_file):
  results = os.path.join(pargenes_dir, "mlsearch_run", "results")
  writer = open(output_file, "w")
  total = 0.0
  for family in sorted(os.listdir(results)):
    logs = os.path.join(results, family, family + ".raxml.log")
    try:
      for line in open(logs).readlines():
        if (line.startswith("Final LogLikelihood:")):
          ll = float(line.split()[2][:-1])
          total += ll
          writer.write(str(ll) + "\n")
    except:
      print("Failed to extract ll for famiy " + family)
  writer.write("\n" + str(total))
  print("Total ll=" + str(total))    

def export_pargenes_trees(pargenes_dir, subst_model, starting_trees, bs_trees, datadir):
  families_dir = os.path.join(datadir, "families")
  
  ok_family = {}
  for family in os.listdir(families_dir):
    ok_family[family] = False
  # tca scores
  concatenated_dir = os.path.join(pargenes_dir, "concatenated_bootstraps")
  fixed_subst_model = subst_model.replace(".", "")
  if (os.path.isdir(concatenated_dir)):
    for concatenation in os.listdir(concatenated_dir):
      family = "_".join(concatenation.split("_")[:-1]) # remove everything after the last
      tca_score = 0.0
      try:
        tca_score = tca.get_tca(os.path.join(concatenated_dir, concatenation), family)
      except:
        print("failed to extract tca score for " + concatenation)
        continue
      try:
        output = os.path.join(families_dir, family, "tca.txt")
        open(output, "w").write(str(tca_score))
      except:
        continue
    # support trees
    support_trees_dir = os.path.join(pargenes_dir, "supports_run", "results")
    for support_tree in os.listdir(support_trees_dir):
      if (not support_tree.endswith("support")):
        continue
      tbe = False
      if ("support.tbe" in support_tree):
        tbe = True
      family = "_".join(support_tree.split("_")[:-1]) # remove everything after the last _
      old_raxml_tree = os.path.join(support_trees_dir, support_tree)
      new_raxml_tree = fam.get_raxml_tree(datadir, fixed_subst_model, family, starting = starting_trees, bstrees = bs_trees, tbe = tbe)
      print(old_raxml_tree)
      print(new_raxml_tree)
      shutil.copyfile(old_raxml_tree, new_raxml_tree)
      ok_family[family] = True

  # ml trees
  ml_trees_dir = os.path.join(pargenes_dir, "mlsearch_run", "results")
  for family in os.listdir(ml_trees_dir):
    trees_file = os.path.join(ml_trees_dir, family, "sorted_ml_trees.newick")
    best_model = os.path.join(ml_trees_dir, family, family + ".raxml.bestModel")
    ml_tree_file = os.path.join(ml_trees_dir, family, family + ".raxml.bestTree")
    if (not os.path.isfile(trees_file)):
      trees_file = ml_tree_file
    family = "_".join(family.split("_")[:-1]) # remove everything after the last _
    new_raxml_trees = fam.get_raxml_multiple_trees(datadir, fixed_subst_model, family, starting_trees)
    new_raxml_tree = fam.get_raxml_tree(datadir, fixed_subst_model, family, starting = starting_trees)
    if (bs_trees == 0):
      shutil.copyfile(ml_tree_file, new_raxml_tree)
      ok_family[family] = True
    try:
      shutil.copyfile(trees_file, new_raxml_trees)
    except:
      print("Cannot copy " + trees_file + " to " + new_raxml_trees)
      pass
    new_best_model = fam.get_raxml_best_model(datadir, fixed_subst_model, family)
    shutil.copyfile(best_model, new_best_model)
  out_likelihood = os.path.join(pargenes_dir, "likelihoods.txt")
  gather_likelihoods(pargenes_dir, out_likelihood)
  print("Likelihoods stored into " + out_likelihood)
  # clean
  garbage_dir = os.path.join(datadir, "garbage")
  try:
    os.makedirs(garbage_dir)
  except:
    pass
  for family in os.listdir(families_dir):
    if (not ok_family[family]):
      print("Cleaning family " + family)
      shutil.move(os.path.join(families_dir, family), garbage_dir)

def get_family_dimensions(datadir, subst_model, pargenes_dir = "pargenes", default_if_fail = False):
  dimensions = {}
  try:
    pargenes_dir = fam.get_run_dir(datadir, subst_model, pargenes_dir)
    parsing_results = os.path.join(pargenes_dir, "parse_run", "results")
    for family in os.listdir(parsing_results):
      log_file = os.path.join(parsing_results, family, family + ".raxml.log")
      unique_sites = 1
      taxa = 1
      lines = open(log_file).readlines()
      for line in lines:
        if "Alignment comprises" in line:
          unique_sites = int(line.split(" ")[5])
        if "taxa and" in line:
          taxa = int(line.split(" ")[4])
      dimensions["_".join(family.split("_")[:-1])] = [unique_sites, taxa]
  except:
    if (default_if_fail):
      print("Setting default dimensions...")
      for family in fam.get_families_list(datadir):
        dimensions[family] = [1, 1]
    else:
      print("from get_family_dimensions...")
      raise
  print("return dimensions")
  return dimensions

  

def run_pargenes_and_extract_trees(datadir, subst_model, rand_trees, parsi_trees, bs_trees, cores, pargenes_dir = "pargenes", extract_trees = True, restart = False, arguments = []):
  saved_metrics_key = "raxml-ng"
  suffix = ""
  if (rand_trees > 0):
    suffix += "-rand" + str(rand_trees)
  if (parsi_trees > 0):
    suffix += "-parsi" + str(parsi_trees)
  if ("--constrain-search" in arguments):
    suffix += "-constrained"
  saved_metrics_key += suffix
  pargenes_dir += suffix
  #if (pargenes_dir != "pargenes"):
  #  saved_metrics_key = pargenes_dir
  pargenes_dir = fam.get_run_dir(datadir, subst_model, pargenes_dir)
  if (not restart):
    shutil.rmtree(pargenes_dir, True)
  start = time.time()
  run_pargenes(datadir, pargenes_dir, subst_model, rand_trees, parsi_trees, bs_trees, cores, arguments)
  saved_metrics.save_metrics(datadir, fam.get_run_name(saved_metrics_key, subst_model), (time.time() - start), "runtimes") 
  lb = fam.get_lb_from_run(os.path.join(pargenes_dir, "mlsearch_run"))
  saved_metrics.save_metrics(datadir, fam.get_run_name(saved_metrics_key, subst_model), (time.time() - start) * lb, "seqtimes") 
  if (extract_trees):
    export_pargenes_trees(pargenes_dir, subst_model, rand_trees + parsi_trees, bs_trees, datadir)

if __name__ == "__main__":
  if (len(sys.argv) < 8):
    print("syntax: python " + os.path.basename(__file__) + " datadir subst_model rand_trees parsi_trees bs_trees cores restart")
    sys.exit(1)
  dataset = sys.argv[1]
  subst_model = sys.argv[2]
  rand_trees= int(sys.argv[3])
  parsi_trees = int(sys.argv[4])
  bs_trees = int(sys.argv[5])
  cores = int(sys.argv[6])
  restart = int(sys.argv[7]) == 1
  arguments = sys.argv[8:]
  run_pargenes_and_extract_trees(dataset, subst_model, rand_trees, parsi_trees, bs_trees, cores, "pargenes", True, restart, arguments)

  
