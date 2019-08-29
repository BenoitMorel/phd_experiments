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

def fix_null_support_values(input_tree, output_tree):
  command = []
  command.append("sed")
  command.append("s/)0:/)1:/g")
  command.append(input_tree)
  with open(output_tree, "w") as writer:
    subprocess.check_call(command, stdout=writer)


def generate_scheduler_commands_file(datadir, subst_model, threshold, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  speciesTree = fam.get_species_tree(datadir)
  family_dimensions = run_pargenes.get_family_dimensions(os.path.abspath(datadir), subst_model)
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_dir = fam.get_family_path(datadir, family)
      mapping_file = fam.get_treerecs_mappings(datadir, family)
      eccetera_dir = fam.get_family_misc_dir(datadir, family)
      eccetera_output = "eccetera." + subst_model + "."
      input_tree = os.path.join(fam.get_family_misc_dir(datadir, family), "eccetera_input." + subst_model + ".geneTree")
      fix_null_support_values(fam.get_raxml_tree(datadir, subst_model, family), input_tree)
      command = []
      command.append(family)
      command.append("1")
      if (family in family_dimensions):
        dim = family_dimensions[family][1] * family_dimensions[family][0]
        command.append(str(dim))
      else:
        command.append("1")
      command.append("species.file=" + speciesTree)
      command.append("gene.file=" + input_tree)
      command.append("output.dir=" + eccetera_dir)
      command.append("output.prefix=" + eccetera_output)
      command.append("resolve.trees=1")
      command.append("print.newick=1")
      command.append("collapse.mode=1")
      command.append("collapse.threshold=" + str(threshold))
      command.append("compute.TD=false")
      command.append("dated=2")
      if(os.path.isfile(mapping_file)):
        command.append("gene.mapping.file=" + mapping_file)
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def extract_eccetera_trees(datadir, subst_model):
  families_dir = os.path.join(datadir, "families")
  for family in os.listdir(families_dir):
    ecceteraTree = os.path.join(families_dir, family, "misc", "eccetera." + subst_model + ".geneTree")
    shutil.copyfile(ecceteraTree, fam.get_eccetera_tree(datadir, subst_model, family)) 

def run_eccetera_on_families(datadir, subst_model, threshold, cores):
  output_dir = fam.get_run_dir(datadir, subst_model, "eccetera_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, threshold, cores, output_dir)
  
  start = time.time()
  exp.run_with_scheduler(exp.eccetera_exec, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, fam.get_run_name("eccetera", subst_model), (time.time() - start), "runtimes") 
  lb = fam.get_lb_from_run(output_dir)
  saved_metrics.save_metrics(datadir, fam.get_run_name("eccetera", subst_model), (time.time() - start) * lb, "seqtimes") 
  extract_eccetera_trees(datadir, subst_model)
  
if (__name__== "__main__"):
  max_args_number = 5
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_eccetera datadir subst_model threshold cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  threshold = int(sys.argv[3])
  cores = int(sys.argv[4])

  run_eccetera_on_families(datadir, subst_model, threshold, cores)


