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


def short_str(num):
  if (num >= 1000000):
    return str(int(num/1000000)) + "M"
  elif (num >= 1000):
    return str(int(num/1000)) + "K"
  else:
    return str(int(num))

class MrbayesInstance():
  def __init__(self, datadir, subst_model, runs = 4, chains = 2, generations = 1000000, frequency = 1000, burnin = 100):
    self.datadir = os.path.abspath(datadir)
    self.subst_model = subst_model
    self.runs = runs
    self.chains = chains
    self.generations = generations
    self.frequency = frequency
    self.burnin = burnin
    basename = self.get_tag("mrbayes_run")
    self.output_dir = fam.get_run_dir(self.datadir, self.subst_model, basename)
    self.output_dir = os.path.abspath(self.output_dir)

  def get_tag(self, prefix = "mrbayes"):
    tag = prefix
    tag += "-r" + short_str(self.runs)
    tag += "-c" + short_str(self.chains)
    tag += "-g" + short_str(self.generations)
    tag += "-f" + short_str(self.frequency)
    tag += "-b" + short_str(self.burnin)
    return tag

  def save_parameters(self):
    output_dir = self.output_dir
    out = os.path.join(output_dir, "mrbayes_parameters.txt")
    with open(out, "w") as writer:
      writer.write("Parameters:")
      writer.write(" " + self.datadir)
      writer.write(" " + str(self.subst_model))
      writer.write(" " + str(self.runs))
      writer.write(" " + str(self.chains))
      writer.write(" " + str(self.generations))
      writer.write(" " + str(self.frequency))
      writer.write(" " + str(self.burnin))
      
  def get_treelist(self, family):
    return os.path.join(self.output_dir, "results", family, family + ".treelist")

  def generate_config_file(self, output_config_file, nexus_alignment, subst_model, seed, output_prefix):
    append = "no"
    ckp = output_prefix + r".ckp~"
    parsi_mode = True
    if (os.path.isfile(ckp)):
      append = "yes"
    with open(output_config_file, "w") as writer:
      writer.write("\tbegin mrbayes;\n")
      writer.write("\tset seed=" + str(seed) + ";\n")
      writer.write("\tset autoclose=yes nowarn=yes;\n")
      writer.write("\texecute " + nexus_alignment + ";\n")
      writer.write(sequence_model.get_mrbayes_preset_line(subst_model))
      writer.write(sequence_model.get_mrbayes_lset_line(subst_model))
      if (parsi_mode):
        writer.write("\tpropset ParsSPR(Tau,V)$prob=0;\n")
        writer.write("\tpropset NNI(Tau,V)$prob=0;\n")
        writer.write("\tpropset   ExtSPR(Tau,V)$prob=0;\n")
        writer.write("\tpropset ExtTBR(Tau,V)$prob=0;\n")
        writer.write("\tpropset ParsSPR1(Tau,V)$prob=12;\n")
        writer.write("\tpropset ParsTBR1(Tau,V)$prob=6;\n")
      writer.write("\tmcmc nruns=1" + " nchains=" + str(self.chains) + " ngen=" + str(self.generations) + " samplefreq=" + str(self.frequency) + " file=" + output_prefix + " append=" + append + ";\n")
      writer.write("end;")

  def remove_mrbayes_run(self):
    output_dir = os.path.abspath(self.output_dir)
    shutil.rmtree(os.path.join(output_dir, "results"), True)

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

def generate_commands_file(instance, cores, prefix_species):
  output_dir = instance.output_dir
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  exp.mkdir(results_dir)
  family_dimensions = {}
  datadir = instance.datadir
  try:
    family_dimensions = run_pargenes.get_family_dimensions(datadir, instance.subst_model)
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
      for run in range(0, instance.runs):
        mrbayes_config = os.path.join(mrbayes_family_dir, "mrbayes_config_run" + str(run) + "." + instance.subst_model + ".nex")
        output_prefix = os.path.join(mrbayes_family_dir, family) + str(run)
        seed = run + 42
        instance.generate_config_file(mrbayes_config, nexus_alignment, instance.subst_model, seed, output_prefix)
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

#def get_reelist(datadir, subst_model, family):
#  return os.path.join(get_mrbayes_output_dir(datadir, subst_model), "results", family, family + ".treelist")


def extract_mrbayes_family(futures_params):
  instance, family = futures_params
  datadir = instance.datadir
  subst_model = instance.subst_model
  burnin = instance.burnin
  tag = instance.get_tag()
  mrbayes_dir = instance.output_dir
  family_misc_dir = fam.get_family_misc_dir(datadir, family) 
  #output = instance.get_treelist(family)
  output = fam.build_gene_tree_path(datadir, subst_model, family, tag)
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

def extract_mrbayes_results(instance):
  start = time.time()
  print("Extracting mrbayes results...")
  futures_params = []
  for family in fam.get_families_list(instance.datadir):
    futures_params.append((instance, family))
  with concurrent.futures.ProcessPoolExecutor() as executor:
    executor.map(extract_mrbayes_family, futures_params)

  print("Finished extracting mrbayes results " + str(time.time() - start) + "s")
  sys.stdout.flush()

def run_mrbayes_on_families(instance, cores, do_continue = False, prefix_species = False):
  cwd = os.getcwd()
  try:
    if (not do_continue):
      shutil.rmtree(instance.output_dir, True)
    try:
      exp.mkdir(instance.output_dir)
    except:
      pass
    print("chdir " + instance.output_dir)
    os.chdir(instance.output_dir)
    instance.output_dir = os.path.relpath(instance.output_dir)
    instance.datadir = os.path.relpath(instance.datadir)
    instance.save_parameters()
    commands = generate_commands_file(instance, cores, prefix_species)
    start = time.time()
    exp.run_with_scheduler(exp.mrbayes_exec, commands, "onecore", cores, instance.output_dir, "logs.txt")   
    tag = instance.get_tag()
    saved_metrics.save_metrics(instance.datadir, fam.get_run_name(tag, instance.subst_model), (time.time() - start), "runtimes") 
    lb = fam.get_lb_from_run(instance.output_dir)
    saved_metrics.save_metrics(instance.datadir, fam.get_run_name(tag, instance.subst_model), (time.time() - start) * lb, "seqtimes") 
    print("Finished running mrbayes after " + str(time.time() - start) + "s")
    sys.stdout.flush()
    extract_mrbayes_results(instance)
  finally:
    os.chdir(cwd)

if (__name__== "__main__"):
  if len(sys.argv) != 10:
    print("Syntax error: python run_mrbayes.py datadir subst_model runs chains generations frequency burnin cores continue{0,1}")
    print(len(sys.argv))
    sys.exit(0)

  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  runs = int(sys.argv[3])
  chains = int(sys.argv[4])
  generations = int(sys.argv[5])
  frequency = int(sys.argv[6])
  burnin = int(sys.argv[7])
  cores = int(sys.argv[8])
  do_continue = int(sys.argv[9]) > 0
  instance = MrbayesInstance(datadir, subst_model, runs, chains, generations, frequency, burnin) 
  run_mrbayes_on_families(instance, cores, do_continue)
  
