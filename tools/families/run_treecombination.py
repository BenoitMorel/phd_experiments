

import sys
import os
import fam
import random
import subprocess
import shutil
import time
import saved_metrics
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
import experiments as exp
import run_raxml_supportvalues

def generate_scheduler_commands_file(datadir, subst_model, tree_number, cores, output_dir, samples):
  results_dir = os.path.join(output_dir, "results")
  samples_dir = os.path.join(output_dir, "samples")
  os.makedirs(results_dir)
  if (samples != tree_number):
    os.makedirs(samples_dir)
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  #family_dimensions = run_raxml_supportvalues.get_family_dimensions(os.path.abspath(datadir), subst_model)
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      prefix = os.path.join(results_dir, family)
      alignment = fam.get_alignment(datadir, family) 
      tree_path = fam.get_raxml_multiple_trees(datadir, "GTR+G", family, starting = tree_number)
      if (tree_number != samples):
        lines = open(tree_path).readlines()
        random.shuffle(lines)
        tree_path = os.path.join(samples_dir, family + ".newick")
        with open(tree_path, "w") as tree_writer:
          for i in range(0, samples):
            tree_writer.write(lines[i].strip())
            if (i != samples - 1):
              tree_writer.write("\n")
      command = []
      command.append(family)
      command.append("1")
      #if (family in family_dimensions):
        #dim = family_dimensions[family][1] * family_dimensions[family][0]
        #command.append(str(dim))
      #else:
      command.append("1")
      command.append(tree_path)
      command.append(alignment)
      command.append(subst_model)
      command.append(prefix)
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def extract_trees(datadir, output_dir, subst_model, tree_number, samples):
  key = "treecombination" + str(tree_number) 
  if (samples != tree_number):
    key = key + "_s" + str(samples)
  results_dir = os.path.join(output_dir, "results")
  families_dir = os.path.join(datadir, "families")
  for family in os.listdir(families_dir):
    prefix = os.path.join(results_dir, family)
    src = prefix + ".newick"
    dest = fam.build_gene_tree_path(datadir, subst_model, family, key)
    shutil.copyfile(src, dest)


def run_treecombination(datadir, subst_model, tree_number, cores, samples):
  key = "treecombination" + str(tree_number) 
  key = key + "_s" + str(samples)
  output_dir = fam.get_run_dir(datadir, subst_model, key)
  do_run = True
  if (do_run):
    shutil.rmtree(output_dir, True)
    os.makedirs(output_dir)
    scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, tree_number, cores, output_dir, samples)
  
    start = time.time()
    exp.run_with_scheduler(exp.treecombine_exec, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
    saved_metrics.save_metrics(datadir, fam.get_run_name(key, subst_model), (time.time() - start), "runtimes") 
    extract_trees(datadir, output_dir, subst_model, tree_number, samples)    


if __name__ == "__main__":
  if (len(sys.argv) < 5):
    print("syntax: python " + os.path.basename(__file__) + " datadir subst_model tree_number cores [samples]")
    sys.exit(1)
  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  tree_number= int(sys.argv[3])
  cores = int(sys.argv[4])
  samples = tree_number
  random.seed(42)
  if (len(sys.argv) > 5):
    samples = int(sys.argv[5])

  dataset = os.path.basename(os.path.normpath(datadir))
  if (dataset == datadir):
    datadir = fam.get_datadir(dataset)
  
  run_treecombination(datadir, subst_model, tree_number, cores, samples)

  

