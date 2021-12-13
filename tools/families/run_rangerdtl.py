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
import sequence_model

def create_ranger_input(species_tree_file, gene_tree_file, output_file):
  with open(output_file, "w") as writer:
    writer.write(open(species_tree_file).read() + "\n")
    writer.write(open(gene_tree_file).read())

def generate_scheduler_commands_file(datadir, subst_model, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  speciesTree = fam.get_species_tree(datadir)
  family_dimensions = run_pargenes.get_family_dimensions(os.path.abspath(datadir), subst_model)
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_dir = fam.get_family_path(datadir, family)
      rangerdtl_dir = fam.get_family_misc_dir(datadir, family)
      gene_tree = fam.get_raxml_tree(datadir, subst_model, family)
      ranger_input = os.path.join(rangerdtl_dir, "input.newick")
      command = []
      command.append(family)
      command.append("1")
      if (family in family_dimensions):
        dim = family_dimensions[family][1] * family_dimensions[family][0]
        command.append(str(dim))
      else:
        command.append("1")
      command.append(exp.rangerdtl_exec)
      command.append("--seed")
      command.append("42")
      command.append("-g")
      command.append(fam.get_raxml_tree(datadir, subst_model, family))
      command.append("-s")
      command.append(speciesTree)
      command.append("-o")
      command.append(rangerdtl_output)
      command.append("-a")
      command.append(alignment_descriptor)
      command.append("-t")
      command.append("all")
      command.append("--ale-evaluation")
      command.append("-T")
      command.append("7")
      command.append("--select-best-tree")
      if(os.path.isfile(mapping_file)):
        command.append("-S")
        command.append(mapping_file)

      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def extract_rangerdtl_trees(datadir, subst_model):
  families_dir = os.path.join(datadir, "families")
  for family in os.listdir(families_dir):
    rangerdtlTree = os.path.join(families_dir, family, "misc", "rangerdtl_output." + subst_model + ".newick.best")
    if (os.path.isfile(rangerdtlTree)):
      lines = open(rangerdtlTree).readlines()
      with open(fam.get_rangerdtl_tree(datadir, subst_model, family), "w") as writer:
        for line in lines:
          if (not line.startswith(">")):
            writer.write(line)
    else:
      print("Warning: no rangerdtl tree for family " + family)  

def run_rangerdtl_on_families(datadir, subst_model, cores):
  output_dir = fam.get_run_dir(datadir, subst_model, "rangerdtl_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, cores, output_dir)
  
  start = time.time()
  exp.run_with_scheduler(exp.rangerdtl_exec, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, fam.get_run_name("rangerdtl", subst_model), (time.time() - start), "runtimes") 
  lb = fam.get_lb_from_run(output_dir)
  saved_metrics.save_metrics(datadir, fam.get_run_name("rangerdtl", subst_model), (time.time() - start) * lb, "seqtimes") 
  extract_rangerdtl_trees(datadir, subst_model)
  
if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_rangerdtl datadir subst_model cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  cores = int(sys.argv[3])

  run_rangerdtl_on_families(datadir, subst_model, cores)


