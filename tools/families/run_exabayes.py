import os
import sys
import subprocess
import shutil
import time
import concurrent.futures
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "scheduler"))
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "msa_edition"))
import scheduler
import saved_metrics
import experiments as exp
import msa_converter
import rename_leaves


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

def generate_config_file(config_file, num_gen, sampling_freq, runs, chains):
  with open(config_file, "w") as writer:
    writer.write("begin run;\n")
    writer.write("  numGen " + str(num_gen) + "\n")
    writer.write("  samplingFreq " + str(sampling_freq) + "\n")
    writer.write("  numCoupledChains " + str(runs) + "\n")
    writer.write("  numRuns " + str(chains) + "\n")
    writer.write("  convergenceCriterion none\n")
    writer.write("end;")

def generate_exabayes_commands_file(dataset_dir, generations, frequency, runs, chains, is_dna, cores, output_dir):
  families_dir = os.path.join(dataset_dir, "families")
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  os.makedirs(results_dir)
  i = 1
  with open(scheduler_commands_file, "w") as writer:
    for family in sorted(os.listdir(families_dir)): # sorted to ensure seeds are unique
# indeed, need different seeds for each exabayes run, because exabayes (old PLL) creates 
# a temporary file whose name comes from rand(). Different seeds avoid collisions
      family_dir = os.path.join(families_dir, family)
      exabayes_family_dir = os.path.join(results_dir, family)
      os.makedirs(exabayes_family_dir)
      exabayes_config = os.path.join(exabayes_family_dir, "exa_config.nex")
      generate_config_file(exabayes_config, generations, frequency, runs, chains)
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
        command.append("PROT")
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



def extract_exabayes_family(params):
  family, exabayes_run_dir, families_dir, burnin = params
  family_misc_dir = os.path.join(families_dir, family, "misc")
  output = os.path.join(family_misc_dir, family + ".treelist")
  with open(output, "w") as writer:
    d = os.path.join(exabayes_run_dir, "results", family)
    for base_file in os.listdir(d):
      if (not "ExaBayes_topologies" in base_file):
        continue
      topologies = os.path.join(d, base_file)
      lines = open(topologies).readlines()
      is_translation = False
      translator = {}
      tree_index = 0
      for line in lines:
        if ("tree gen" in line):
          tree_index += 1
          if (tree_index <= burnin):
            continue
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

def extract_exabayes_results(exabayes_run_dir, families_dir, burnin):
  start = time.time()
  print("Extracting exabayes results...")
  params = []
  for family in os.listdir(families_dir):
    params.append((family, exabayes_run_dir, families_dir, burnin))

  with concurrent.futures.ProcessPoolExecutor() as executor:
    executor.map(extract_exabayes_family, params)
  shutil.rmtree( os.path.join(exabayes_run_dir, "results"))
  print("Finished extracting exabayes results " + str(time.time() - start) + "s")
  sys.stdout.flush()
  
def get_exabayes_output_dir(dataset_dir):
  return os.path.join(dataset_dir, "runs", "exabayes_run")
  

def run_exabayes_on_families(dataset_dir, generations, frequency, runs, chains, burnin, is_dna, cores):
  fam.init_dataset_dir(dataset_dir)
  output_dir = get_exabayes_output_dir(dataset_dir)
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_exabayes_commands_file(dataset_dir, generations, frequency, runs, chains, is_dna, cores, output_dir)
  start = time.time()
  scheduler.run_scheduler(scheduler_commands_file, exp.exabayes_exec, cores, output_dir, "exabayes_run.logs")
  saved_metrics.save_metrics(dataset_dir, "ExaBayes", (time.time() - start), "runtimes") 
  print("Finished running exabayes after " + str(time.time() - start) + "s")
  sys.stdout.flush()
  extract_exabayes_results(output_dir, os.path.join(dataset_dir, "families"), burnin)

def clean_exabayes(dataset_dir):
  families_dir = os.path.join(dataset_dir, "families")
  for family in os.listdir(families_dir):
    family_dir = os.path.join(families_dir, family)
    family_misc_dir = os.path.join(family_dir, "misc")
    treelist = os.path.join(family_misc_dir, family + ".treelist")
    os.remove(os.path.join(treelist))
  shutil.rmtree(get_exabayes_output_dir(dataset_dir))

if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_exabayes.py dataset_dir is_dna cores.")
    sys.exit(0)

  dataset_dir = sys.argv[1]
  generations = int(sys.argv[2])
  frequency = int(sys.argv[3])
  runs = int(sys.argv[4])
  chains = int(sys.argv[5])
  cores = int(sys.argv[6])
  burnin = int(sys.argv[7])
  is_dna = int(sys.argv[8]) != 0
  run_exabayes_on_families(dataset_dir, generations, frequency, runs, chains, burnin, is_dna, cores)


#
