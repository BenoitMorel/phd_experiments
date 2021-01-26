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
import run_raxml_supportvalues as run_pargenes

def get_mrbayes_output_dir(datadir, subst_model):
  return fam.get_run_dir(datadir, subst_model, "mrbayes_run")

def get_treelist(datadir, subst_model, family):
    return os.path.join(get_mrbayes_output_dir(datadir, subst_model), "results", family, family + ".treelist")


def generate_config_file(mrbayes_config, generations, frequency, chains, nexus_alignment, subst_model, seed, output_prefix):
  with open(mrbayes_config, "w") as writer:
    writer.write("\tbegin mrbayes;\n")
    writer.write("\tset seed=" + str(seed) + ";\n")
    writer.write("\tset autoclose=yes nowarn=yes;\n")
    writer.write("\texecute " + nexus_alignment + ";\n")
    writer.write(sequence_model.get_mrbayes_preset_line(subst_model))
    writer.write(sequence_model.get_mrbayes_lset_line(subst_model))
    writer.write("\tmcmc nruns=1" + " nchains=" + str(chains) + " ngen=" + str(generations) + " samplefreq=" + str(frequency) + " file=" + output_prefix + ";\n")
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

def generate_mrbayes_commands_file(datadir, generations, frequency, runs, chains, subst_model, cores, output_dir, prefix_species):
  output_dir = get_mrbayes_output_dir(datadir, subst_model)
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  exp.mkdir(results_dir)
  family_dimensions = {}
  try:
    family_dimensions = run_pargenes.get_family_dimensions(os.path.abspath(datadir), subst_model)
  except:
    pass
  with open(scheduler_commands_file, "w") as writer:
    for family in sorted(fam.get_families_list(datadir)):
      family_dir = fam.get_family_path(datadir, family)
      mrbayes_family_dir = os.path.join(results_dir, family)
      exp.mkdir(mrbayes_family_dir)
      nexus_alignment = os.path.join(family_dir, "species_prefixed_alignment.nex")
      fasta_alignment = os.path.join(family_dir, "alignment.msa")
      mapping_dictionnary = None
      if (prefix_species):
        mapping_dictionnary = get_mapping_dictionnary(fam.get_mappings(datadir, family))
      msa_converter.msa_convert(fasta_alignment, nexus_alignment, "fasta", "nexus", mapping_dictionnary)
      for run in range(0, runs):
        mrbayes_config = os.path.join(mrbayes_family_dir, "mrbayes_config_run" + str(run) + "." + subst_model + ".nex")
        output_prefix = os.path.join(mrbayes_family_dir, family) + str(run)
        generate_config_file(mrbayes_config, generations, frequency, chains, nexus_alignment, subst_model, run + 42, output_prefix)
        command = []
        command.append(family + "__" + str(run))
        command.append("1")
        if (family in family_dimensions):
          dim = family_dimensions[family][1] * family_dimensions[family][0]
          command.append(str(dim))
        else:
          command.append("1")
        command.append(mrbayes_config)
        writer.write(" ".join(command) + "\n")
  return scheduler_commands_file

def get_reelist(datadir, subst_model, family):
  return os.path.join(get_mrbayes_output_dir(datadir, subst_model), "results", family, family + ".treelist")

def remove_mrbayes_run(datadir, subst_model):
  output_dir = os.path.abspath(get_mrbayes_output_dir(datadir, subst_model))
  shutil.rmtree(os.path.join(output_dir, "results"), True)

def extract_mrbayes_family(params):
  datadir, subst_model, family, mrbayes_dir, burnin = params
  family_mist_dir = fam.get_family_misc_dir(datadir, family) 
  output = get_treelist(datadir, subst_model, family)
  gene_tree_path = fam.build_gene_tree_path(datadir, subst_model, family, "mrbayes")
  output = gene_tree_path
  with open(output, "w") as writer:
    d = os.path.join(mrbayes_dir, "results", family)
    for topologies in os.listdir(d):
      if (not topologies.endswith(".t")):
        continue
      topologies = os.path.join(d, topologies)
      
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
          split = line.split(" ")
          left = split[-2]
          right = split[-1][:-2]
          translator[left] = right

def extract_mrbayes_results(datadir, subst_model, burnin):
  output_dir = get_mrbayes_output_dir(datadir, subst_model)
  start = time.time()
  print("Extracting mrbayes results...")
  params = []
  for family in fam.get_families_list(datadir):
    params.append((datadir, subst_model, family, output_dir, burnin))
  with concurrent.futures.ProcessPoolExecutor() as executor:
    executor.map(extract_mrbayes_family, params)

  print("Finished extracting mrbayes results " + str(time.time() - start) + "s")
  sys.stdout.flush()

def run_mrbayes_on_families(datadir, generations, frequency, runs, chains, burnin, subst_model, cores, prefix_species = False):
  output_dir = os.path.abspath(get_mrbayes_output_dir(datadir, subst_model))
  datadir = os.path.abspath(datadir)
  cwd = os.getcwd()
  try:
    shutil.rmtree(output_dir, True)
    exp.mkdir(output_dir)
    print("chdir " + output_dir)
    os.chdir(output_dir)
    output_dir = os.path.relpath(output_dir)
    datadir = os.path.relpath(datadir)
    parameters = os.path.join(output_dir, "..", "mrbayes_parameters." + subst_model + ".txt")
    open(parameters, "w").write("Parameters: " + datadir + " " + str(generations) + " " + str(frequency) + " " + str(runs) + " " + str(chains) + " " + str(burnin) + " " + str(subst_model) + " " + str(cores))
    scheduler_commands_file = generate_mrbayes_commands_file(datadir, generations, frequency, runs, chains, subst_model, cores, output_dir, prefix_species)
    start = time.time()
    exp.run_with_scheduler(exp.mrbayes_exec, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
    saved_metrics.save_metrics(datadir, fam.get_run_name("mrbayes", subst_model), (time.time() - start), "runtimes") 
    lb = fam.get_lb_from_run(output_dir)
    saved_metrics.save_metrics(datadir, fam.get_run_name("mrbayes", subst_model), (time.time() - start) * lb, "seqtimes") 
    print("Finished running mrbayes after " + str(time.time() - start) + "s")
    sys.stdout.flush()
    extract_mrbayes_results(datadir, subst_model, burnin)
  finally:
    os.chdir(cwd)

if (__name__== "__main__"):
  if len(sys.argv) != 9:
    print("Syntax error: python run_mrbayes.py datadir generations frequency runs chains burnin subst_model cores")
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

