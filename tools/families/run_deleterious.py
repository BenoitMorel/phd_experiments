import os
import sys
import subprocess
import shutil
import time
import saved_metrics
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
import experiments as exp
import fam
import run_raxml_supportvalues as run_pargenes
import ete3
import sequence_model

def add_species_bl(input_species_tree, output_species_tree):
  print(input_species_tree)
  tree = ete3.Tree(input_species_tree, format = 1)
  if (tree.dist == 0.0):
    tree.dist = 1.0
  with open(output_species_tree, "w") as writer:
    writer.write(tree.write(format_root_node = True))

def generate_scheduler_commands_file(datadir, subst_model, cores, output_dir):
  generations = 10000
  thinning = 10
  results_dir = os.path.join(output_dir, "results")
  deleterious_species_tree = os.path.join(output_dir, "speciesTree.newick")
  add_species_bl(fam.get_species_tree(datadir), deleterious_species_tree)
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  family_dimensions = run_pargenes.get_family_dimensions(os.path.abspath(datadir), subst_model, default_if_fail = True)
  print("end of getting family dimensions")
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_dir = fam.get_family_path(datadir, family)
      mapping_file = fam.get_treerecs_mappings(datadir, family)
      results = os.path.join(results_dir, family)
      exp.mkdir(results)
      command = []
      command.append(family)
      command.append("1")
      if (family in family_dimensions):
        dim = family_dimensions[family][1] * family_dimensions[family][0]
        command.append(str(dim))
      else:
        command.append("1")
      command.append("-Xms512m")
      command.append("-Xmx1024m")
      command.append("-jar")
      command.append(exp.jprime_jar)
      command.append("Deleterious")
      command.append(deleterious_species_tree)
      command.append(fam.get_alignment(datadir, family))
      command.append(mapping_file)
      command.append("-i")
      command.append(str(generations))
      command.append("-t")
      command.append(str(thinning))
      command.append("-s")
      command.append("42")
      command.append("--siteratecategories")
      command.append(str(sequence_model.get_gamma_rates(subst_model)))
      command.append("-sm")
      command.append(sequence_model.get_deleterious_model(subst_model))
      command.append("-o")
      command.append(os.path.join(results, "deleterious"))
      writer.write(" ".join(command) + "\n")
  print("wrote command file")
  return scheduler_commands_file
    
def extract_tree_from_run(datadir, subst_model, family):
  output_dir = fam.get_run_dir(datadir, subst_model, "deleterious_run")
  tree_info = os.path.join(output_dir, "results", family, "deleterious.info")
  output_tree = fam.get_deleterious_tree(datadir, subst_model, family)
  try:
    with open(tree_info) as reader:
      while (True):
        line = reader.readline()
        if (not line):
          break
        if (line.endswith(";\n") and not "Host" in line):
          with open(output_tree, "w") as writer:
            writer.write(line.split()[-1])
            return
  except:
    pass
  print("No deleterious tree for family " + family)
  shutil.copyfile(fam.get_raxml_tree(datadir, subst_model, family), output_tree)

def extract_deleterious_trees(datadir, subst_model):
  families_dir = os.path.join(datadir, "families")
  for family in os.listdir(families_dir):
    extract_tree_from_run(datadir, subst_model, family)

def run_deleterious_on_families(datadir, subst_model,  cores):
  output_dir = fam.get_run_dir(datadir, subst_model, "deleterious_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, cores, output_dir)
  
  start = time.time()
  print("before run")
  exp.run_with_scheduler("java", scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  print("after run")
  saved_metrics.save_metrics(datadir, fam.get_run_name("deleterious", subst_model), (time.time() - start), "runtimes") 
  lb = fam.get_lb_from_run(output_dir)
  saved_metrics.save_metrics(datadir, fam.get_run_name("deleterious", subst_model), (time.time() - start) * lb, "seqtimes") 
  extract_deleterious_trees(datadir, subst_model)
  
if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_deleterious datadir subst_model cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  cores = int(sys.argv[3])

  run_deleterious_on_families(datadir, subst_model,  cores)
  extract_deleterious_trees(datadir, subst_model)


