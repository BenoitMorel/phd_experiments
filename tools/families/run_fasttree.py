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
import sequence_model
import ete3

def generate_scheduler_commands_file(datadir, subst_model, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  gamma = False
  base_model = subst_model
  if ("+G" in subst_model):
    gamma = True
    base_model = subst_model.replace("+G", "")
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_dir = fam.get_family_path(datadir, family)
      fasttree_dir = fam.get_family_misc_dir(datadir, family)
      try:
        os.mkdir(fasttree_dir)
      except:
        pass
      fasttree_output = os.path.join(fasttree_dir, "fasttree_output." + subst_model + ".newick")
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      #command.append(exp.fasttree_exec)
      #../FastTree/FastTree -out out.newick -nt -gtr
      if ("WAG" in subst_model):
        command.append("-wag")
      elif ("LG" in subst_model):
        command.append("-lg")
      else:
        command.append("-nt")
      if (gamma):
        command.append("-gamma")
      command.append("-" + base_model.lower())
      command.append("-out")
      command.append(fasttree_output)
      command.append(fam.get_alignment(datadir, family))
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def extract_fasttree_trees(datadir, subst_model):
  families_dir = os.path.join(datadir, "families")
  for family in os.listdir(families_dir):
    fasttreeTree = os.path.join(families_dir, family, "misc", "fasttree_output." + subst_model + ".newick")
    if (os.path.isfile(fasttreeTree)):
      lines = open(fasttreeTree).readlines()
      writerpoly = open(fam.get_fasttreepoly_tree(datadir, subst_model, family), "w")
      with open(fam.get_fasttree_tree(datadir, subst_model, family), "w") as writer:
        for line in lines:
          if (not line.startswith(">")):
            writerpoly.write(line)
            tree = ete3.Tree(line, format = 1)
            tree.resolve_polytomy()
            writer.write(tree.write())
    else:
      print("Warning: no fasttree tree for family " + family)  

def run_fasttree_on_families(datadir, subst_model, cores):
  output_dir = fam.get_run_dir(datadir, subst_model, "fasttree_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, cores, output_dir)
  
  start = time.time()
  exp.run_with_scheduler(exp.fasttree_exec, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, fam.get_run_name("fasttree", subst_model), (time.time() - start), "runtimes") 
  lb = fam.get_lb_from_run(output_dir)
  saved_metrics.save_metrics(datadir, fam.get_run_name("fasttree", subst_model), (time.time() - start) * lb, "seqtimes") 
  extract_fasttree_trees(datadir, subst_model)
  
if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python " + os.path.basename(__file__) + " datadir subst_model cores.")
    sys.exit(0)


  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  cores = int(sys.argv[3])

  run_fasttree_on_families(datadir, subst_model, cores)


