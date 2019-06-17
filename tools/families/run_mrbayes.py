import os
import sys
import subprocess
import shutil
import time
import concurrent.futures
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "msa_edition"))
import saved_metrics
import experiments as exp
import msa_converter
import rename_leaves
import sequence_model

def get_mrbayes_output_dir(datadir, subst_model):
  return fam.get_run_dir(datadir, subst_model, "mrbayes_run")
      

def generate_config_file(mrbayes_config, generations, frequency, runs, chains, nexus_alignment, subst_model, output_prefix):
  with open(mrbayes_config, "w") as writer:
    writer.write("\tbegin mrbayes;\n")
    writer.write("\tset autoclose=yes nowarn=yes;\n")
    writer.write("\texecute " + nexus_alignment + ";\n")
    writer.write(sequence_model.get_mrbayes_preset_line(subst_model))
    writer.write(sequence_model.get_mrbayes_lset_line(subst_model))
    writer.write("\tmcmc nruns=" + str(runs) + " ngen=" + str(generations) + " samplefreq=" + str(frequency) + " file=" + output_prefix + ";\n")
    writer.write("end;")

def get_mapping_dictionnary(mapping_file):
  res = {}
  lines = open(mapping_file).readlines()
  for line in lines:
    split = line[:-1].split(":")
    species = split[0]
    genes = split[1].split(";")
    for gene in genes:
      res[gene] = species
  return res

def generate_mrbayes_commands_file(datadir, generations, frequency, runs, chains, subst_model, cores, output_dir):
  output_dir = get_mrbayes_output_dir(datadir, subst_model)
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  exp.mkdir(results_dir)
  with open(scheduler_commands_file, "w") as writer:
    for family in sorted(fam.get_families_list(datadir)):
      family_dir = fam.get_family_path(datadir, family)
      mrbayes_family_dir = os.path.join(results_dir, family)
      exp.mkdir(mrbayes_family_dir)
      mrbayes_config = os.path.join(mrbayes_family_dir, "mrbayes_config.nex")
      nexus_alignment = os.path.join(family_dir, "species_prefixed_alignment.nex")
      fasta_alignment = os.path.join(family_dir, "alignment.msa")
      mapping_dictionnary = get_mapping_dictionnary(fam.get_mappings(datadir, family))
      msa_converter.msa_convert(fasta_alignment, nexus_alignment, "fasta", "nexus", mapping_dictionnary)
      output_prefix = os.path.join(mrbayes_family_dir, family)
      generate_config_file(mrbayes_config, generations, frequency, runs, chains, nexus_alignment, subst_model, output_prefix)
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append(mrbayes_config)
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
  
def extract_mrbayes_results(datadir, burnin):
  output_dir = get_mrbayes_output_dir(datadir, subst_model)
  print("TODO IMPLEMENT")

def run_mrbayes_on_families(datadir, generations, frequency, runs, chains, burnin, subst_model, cores):
  output_dir = get_mrbayes_output_dir(datadir, subst_model)
  shutil.rmtree(output_dir, True)
  exp.mkdir(output_dir)
  parameters = os.path.join(output_dir, "parameters.txt")
  open(parameters, "w").write("Parameters: " + datadir + " " + str(generations) + " " + str(frequency) + " " + str(runs) + " " + str(chains) + " " + str(burnin) + " " + str(subst_model) + " " + str(cores))
  scheduler_commands_file = generate_mrbayes_commands_file(datadir, generations, frequency, runs, chains, subst_model, cores, output_dir)
  start = time.time()
  exp.run_with_scheduler(exp.mrbayes_exec, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, "mrbayes", (time.time() - start), "runtimes") 
  print("Finished running mrbayes after " + str(time.time() - start) + "s")
  sys.stdout.flush()
  extract_mrbayes_results(datadir, burnin)
  print("TODO SET CHAINS")

if (__name__== "__main__"):
  if len(sys.argv) != 9:
    print("Syntax error: python run_mrbayes.py datadir generations frequency runs chains cores burnin subst_model")
    print(len(sys.argv))
    sys.exit(0)

  datadir = sys.argv[1]
  generations = int(sys.argv[2])
  frequency = int(sys.argv[3])
  runs = int(sys.argv[4])
  chains = int(sys.argv[5])
  burnin = int(sys.argv[6])
  subst_model = sys.argv[7]
  cores = int(sys.argv[8])
  run_mrbayes_on_families(datadir, generations, frequency, runs, chains, burnin, subst_model, cores)


#

