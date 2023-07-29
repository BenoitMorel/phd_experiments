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




def generate_fastme_scheduler_commands_file(datadir, subst_model, is_dna, cores, output_dir):
  results_dir = os.path.join(output_dir, "resultsfastme")
  scheduler_commands_file = os.path.join(output_dir, "fastme_commands.txt")
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
      command.append("-c")
      command.append("--seed")
      command.append("40")
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file

def remove_blank_lines(matrix_file):
  lines = open(matrix_file).readlines()
  with open(matrix_file, "w") as writer:
    for line in lines:
      if line in ['\n', '\r\n']:
        continue
      writer.write(line)

def generate_noah_scheduler_commands_file(datadir, subst_model, is_dna, samples, perturbation, cores, output_dir):
  results_dir = os.path.join(output_dir, "resultsnoah")
  scheduler_commands_file = os.path.join(output_dir, "noah_commands.txt")
  gamma = False
  sp = subst_model.split("+")
  fastme_model = subst_model
  if (len(sp) > 1 and sp[1] == "G"):
    fastme_model = sp[0]
    gamma = True
  run_name = "noah_s" + str(samples) + "_p" + str(perturbation).replace(".", "_")
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      fastme_matrix = fam.get_fastme_distances(datadir, family, subst_model)
      remove_blank_lines(fastme_matrix)
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append("--sequence-file")
      command.append(fastme_matrix)
      command.append("-o")
      tree = fam.build_gene_tree_path(datadir, subst_model, family, run_name)
      command.append(tree)
      command.append("-t")
      command.append(str(samples))
      command.append("--seed")
      command.append("40")
      command.append("-n")
      command.append(str(perturbation))
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file


def run_noah_on_families(datadir, subst_model, is_dna, samples, perturbation, cores):
  run_name = "noah_s" + str(samples) + "_p" + str(perturbation).replace(".", "_")
  output_dir = fam.get_run_dir(datadir, subst_model, run_name + "_run")
  fastme_output_dir = fam.get_run_dir(output_dir, "fastme_scheduler")
  noah_output_dir = fam.get_run_dir(output_dir, "noah_scheduler")

  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  os.makedirs(noah_output_dir)
  os.makedirs(fastme_output_dir)
  fastme_commands_file = generate_fastme_scheduler_commands_file(datadir, subst_model, is_dna, cores, output_dir)
  start = time.time()
  exp.run_with_scheduler(exp.fastme_exec, fastme_commands_file, "onecore", cores, fastme_output_dir, "logs.txt")   
  noah_commands_file = generate_noah_scheduler_commands_file(datadir, subst_model, is_dna, samples, perturbation, cores, output_dir)
  exp.run_with_scheduler(exp.noah_exec, noah_commands_file, "onecore", cores, noah_output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), (time.time() - start), "runtimes") 
  #lb = fam.get_lb_from_run(output_dir)
  #saved_metrics.save_metrics(datadir, fam.get_run_name("noah", subst_model), (time.time() - start) * lb, "seqtimes") 
 

if (__name__== "__main__"):
  max_args_number = 7
  if len(sys.argv) < max_args_number:
    print("Syntax error: python " + os.path.basename(__file__) + "  datadir subst_model is_dna samples perturbation cores.")
    sys.exit(0)
  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  is_dna = sys.argv[3] != "0"
  samples = int(sys.argv[4])
  perturbation = float(sys.argv[5])
  cores = int(sys.argv[6])
  run_noah_on_families(datadir, subst_model, is_dna, samples, perturbation, cores)



