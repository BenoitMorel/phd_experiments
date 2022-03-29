import os
import sys
import subprocess
import shutil
import time
import saved_metrics
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "msa_edition"))
import experiments as exp
import fam
import ete3
import msa_converter




def generate_scheduler_commands_file(datadir, subst_model, is_dna, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  gamma = False
  sp = subst_model.split("+")
  fastme_model = subst_model
  if (len(sp) > 1 and sp[1] == "G"):
    fastme_model = sp[0]
    gamma = True
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_dir = fam.get_family_path(datadir, family)
      fastme_dir = fam.get_family_misc_dir(datadir, family)
      try:
        os.mkdir(fastme_dir)
      except:
        pass
      fastme_output = os.path.join(fastme_dir, "fastme_output." + subst_model + ".newick")
      fastme_matrix = fam.get_fastme_distances(datadir, family, subst_model)
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append("-i")
      phylip = fam.get_alignment_phylip(datadir, family)
      if (not os.path.isfile(phylip)):
        ali = fam.get_alignment(datadir, family)
        msa_converter.msa_convert(ali, phylip, None, "iphylip_relaxed")
      command.append(phylip)
      
      if (is_dna):
        command.append("-d" + fastme_model)
      else:
        command.append("-p" + fastme_model)
      if (gamma):
        command.append("1.0")
      command.append("-o")
      command.append(fastme_output)
      command.append("-O")
      command.append(fastme_matrix)
      command.append("--seed")
      command.append("40")
      command.append("--spr")
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file


def extract_fastme_trees(datadir, subst_model):
  families_dir = os.path.join(datadir, "families")
  valid = 0
  invalid = 0
  for family in os.listdir(families_dir):
    fastmetree = os.path.join(families_dir, family, "misc", "fastme_output." + subst_model + ".newick")
    tree = fam.build_gene_tree_path(datadir, subst_model, family, "fastme")
    fastme_matrix = fam.get_fastme_distances(datadir, family, subst_model)
    if (os.path.isfile(fastmetree) and os.stat(fastmetree).st_size > 0):
      valid += 1
      shutil.copyfile(fastmetree, tree)
    else:
      invalid += 1
      try:
        os.remove(tree)
        os.remove(fastme_matrix)
      except:
        pass
    os.remove(fastmetree) 
  print("Extracted " + str(valid) + " trees")
  if (invalid > 0):
    print("WARNING! " + str(invalid) + " trees were skipped")

def run_fastme_on_families(datadir, subst_model, is_dna, cores):
  output_dir = fam.get_run_dir(datadir, subst_model, "fastme_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, is_dna, cores, output_dir)
  start = time.time()
  exp.run_with_scheduler(exp.fastme_exec, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, fam.get_run_name("fastme", subst_model), (time.time() - start), "runtimes") 
  lb = fam.get_lb_from_run(output_dir)
  saved_metrics.save_metrics(datadir, fam.get_run_name("fastme", subst_model), (time.time() - start) * lb, "seqtimes") 
  extract_fastme_trees(datadir, subst_model)
 

if (__name__== "__main__"):
  max_args_number = 5
  if len(sys.argv) < max_args_number:
    print("Syntax error: python " + os.path.basename(__file__) + "  datadir subst_model is_dna cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)
  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  is_dna = sys.argv[3] != "0"
  cores = int(sys.argv[4])
  run_fastme_on_families(datadir, subst_model, is_dna, cores)


