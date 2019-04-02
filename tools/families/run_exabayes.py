import os
import sys
import subprocess
import shutil
import time
import utils
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "msa_edition"))
import saved_metrics
import experiments as exp
import msa_converter
import rename_leaves

EXA_GEN = 40000
EXA_FREQ = 200

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

def generate_config_file(config_file, num_gen, sampling_freq):
  with open(config_file, "w") as writer:
    writer.write("begin run;\n")
    writer.write("  numGen " + str(num_gen) + "\n")
    writer.write("  samplingFreq " + str(sampling_freq) + "\n")
    writer.write("end;")

def generate_exabayes_commands_file(dataset_dir, is_dna, cores, output_dir):
  families_dir = os.path.join(dataset_dir, "families")
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  os.makedirs(results_dir)
  i = 1
  with open(scheduler_commands_file, "w") as writer:
    for family in os.listdir(families_dir):
      family_dir = os.path.join(families_dir, family)
      exabayes_family_dir = os.path.join(results_dir, family)
      os.makedirs(exabayes_family_dir)
      exabayes_config = os.path.join(exabayes_family_dir, "exa_config.nex")
      generate_config_file(exabayes_config, EXA_GEN, EXA_FREQ)
      phy_alignment = os.path.join(family_dir, "species_prefixed_alignment.phy")
      fasta_alignment = os.path.join(family_dir, "alignment.msa")
      mapping_dictionnary = get_mapping_dictionnary(os.path.join(family_dir, "mapping.link"))
      msa_converter.msa_convert(fasta_alignment, phy_alignment, "fasta", "iphylip_relaxed", mapping_dictionnary)
      run_id = family
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append("-f")
      command.append(phy_alignment)
      command.append("-m")
      if (is_dna):
        command.append("DNA")
      else:
        command.append("LG")
      command.append("-s")
      command.append(str(i))
      command.append("-n")
      command.append(run_id)
      command.append("-c")
      command.append(exabayes_config)
      command.append("-w")
      command.append(exabayes_family_dir)
      writer.write(" ".join(command) + "\n")
      i += 1
  return scheduler_commands_file


def extract_exabayes_results(exabayes_run_dir, families_dir):
  for family in os.listdir(families_dir):
    family_misc_dir = os.path.join(families_dir, family, "misc")
    try:
      os.makedirs(family_misc_dir)
    except:
      pass
    topologies = os.path.join(exabayes_run_dir, "results", family, "ExaBayes_topologies.run-0." + family)
    lines = open(topologies).readlines()
    output = os.path.join(family_misc_dir, family + ".treelist")
    is_translation = False
    translator = {}
    with open(output, "w") as writer:
      for line in lines:
        if ("tree gen" in line):
          is_translation = False
          tree = line.split(" ")[-1]
          writer.write(rename_leaves.rename_leaves(tree, translator) + "\n")
          continue
        if ("translate" in line):
          is_translation = True
          continue
        if (is_translation):
          split = line.split("\t")
          left = split[1]
          right = split[2][:-2]
          translator[left] = right

def run_exabayes_on_families(dataset_dir, is_dna, cores):
  output_dir = os.path.join(dataset_dir, "exabayes_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_exabayes_commands_file(dataset_dir, is_dna, cores, output_dir)
  start = time.time()
  utils.run_scheduler(scheduler_commands_file, exp.exabayes_exec, cores, output_dir, "exabayes_run.logs")
  saved_metrics.save_metrics(dataset_dir, "ExaBayes", (time.time() - start), "runtimes") 
  extract_exabayes_results(output_dir, os.path.join(dataset_dir, "families"))

if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_exabayes.py dataset_dir is_dna cores.")
    sys.exit(0)

  dataset_dir = sys.argv[1]
  is_dna = int(sys.argv[2]) != 0
  cores = int(sys.argv[3])
  run_exabayes_on_families(dataset_dir, is_dna, cores)


#
