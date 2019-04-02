import os
import sys
import subprocess
import shutil
import time
import utils
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "msa_edition"))
import saved_metrics
import experiments as exp
import msa_converter

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

def generate_phylobayes_commands_file(dataset_dir, cores, output_dir):
  families_dir = os.path.join(dataset_dir, "families")
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  os.makedirs(results_dir)
  with open(scheduler_commands_file, "w") as writer:
    for family in os.listdir(families_dir):
      family_dir = os.path.join(families_dir, family)
      phylobayes_dir = os.path.join(family_dir, "phylobayes")
      phy_alignment = os.path.join(os.path.join(family_dir, "species_prefixed_alignment.phy"))
      fasta_alignment = os.path.join(os.path.join(family_dir, "alignment.msa"))
      mapping_dictionnary = get_mapping_dictionnary(os.path.join(family_dir, "mapping.link"))
      if (not os.path.isfile(phy_alignment)):
        msa_converter.msa_convert(fasta_alignment, phy_alignment, "fasta", "iphylip_relaxed", mapping_dictionnary)
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append("-d")
      command.append(phy_alignment)
      command.append("-x")
      command.append("1")
      command.append("4000")
      command.append(os.path.join(results_dir, family))
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file


def extract_phylobayes_results(phylobayes_run_dir, families_dir):
  for family in os.listdir(families_dir):
    family_misc_dir = os.path.join(families_dir, family, "misc")
    try:
      os.makedirs(family_misc_dir)
    except:
      pass
    treelist = os.path.join(phylobayes_run_dir, "results", family + ".treelist")
    shutil.copyfile(treelist, os.path.join(family_misc_dir, family + ".treelist"))

def run_phylobayes_on_families(dataset_dir, cores):
  output_dir = os.path.join(dataset_dir, "phylobayes_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_phylobayes_commands_file(dataset_dir, cores, output_dir)
  start = time.time()
  utils.run_scheduler(scheduler_commands_file, exp.phylobayes_exec, cores, output_dir, "phylobayes_run.logs")
  saved_metrics.save_metrics(dataset_dir, "PhyloBayes", (time.time() - start), "runtimes") 
  extract_phylobayes_results(output_dir, os.path.join(dataset_dir, "families"))

if (__name__== "__main__"):
  max_args_number = 3
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_phylobayes.py dataset_dir cores.")
    sys.exit(0)

  dataset_dir = sys.argv[1]
  cores = int(sys.argv[2])
  run_phylobayes_on_families(dataset_dir, cores)


#
