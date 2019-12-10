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
import ete3

def root_tree(input_tree, output_tree):
  tree = ete3.Tree(input_tree, 1)
  tree.resolve_polytomy()
  open(output_tree, "w").write(tree.write())

def spaces_to_tabs(input_file, output_file):
  f = open(input_file).read()
  open(output_file, "w").write(f.replace(" ", "\t"))

def generate_scheduler_commands_file(datadir, subst_model, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  speciesTree = fam.get_species_tree(datadir)
  family_dimensions = run_pargenes.get_family_dimensions(os.path.abspath(datadir), subst_model)
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_dir = fam.get_family_path(datadir, family)
      treefix_dir = fam.get_family_misc_dir(datadir, family)
      treefix_startingtree = os.path.join(treefix_dir, family + ".treefix.newick.in")
      treefix_alignment = os.path.abspath(os.path.join(treefix_dir, family + ".treefix.msa"))
      treefix_mappping = os.path.join(treefix_dir, family + ".treefix.map")
      
      #exp.relative_symlink(fam.get_alignment(datadir, family), treefix_alignment)
      shutil.copyfile(fam.get_alignment(datadir, family), treefix_alignment)
      root_tree(fam.get_raxml_tree(datadir, subst_model, family), treefix_startingtree)
      spaces_to_tabs(fam.get_treerecs_mappings(datadir, family), treefix_mappping)
      command = []
      command.append(family)
      command.append("1")
      if (family in family_dimensions):
        dim = family_dimensions[family][1] * family_dimensions[family][0]
        command.append(str(dim))
      else:
        command.append("1")
      #command.append(exp.treefix_exec)
      command.append("-s")
      command.append(speciesTree)
      command.append("-S")
      command.append(treefix_mappping)
      command.append("-A")
      command.append(".treefix.msa")
      print(treefix_alignment)
      command.append("-o")
      command.append(".treefix.newick.in")
      command.append("-n")
      command.append(".treefix.newick.out")
      command.append(treefix_startingtree)
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def extract_treefix_trees(datadir, subst_model):
  families_dir = os.path.join(datadir, "families")
  for family in os.listdir(families_dir):
    treefix_dir = fam.get_family_misc_dir(datadir, family)
    treefix_output_tree = os.path.join(treefix_dir, family + ".treefix.newick.out")
    shutil.copy(treefix_output_tree, fam.get_treefix_tree(datadir, subst_model, family))

def run_treefix_on_families(datadir, subst_model, cores):
  output_dir = fam.get_run_dir(datadir, subst_model, "treefix_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, cores, output_dir)
  start = time.time()
  exp.run_with_scheduler(exp.treefix_exec, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, fam.get_run_name("treefix", subst_model), (time.time() - start), "runtimes") 
  lb = fam.get_lb_from_run(output_dir)
  saved_metrics.save_metrics(datadir, fam.get_run_name("treefix", subst_model), (time.time() - start) * lb, "seqtimes") 
  extract_treefix_trees(datadir, subst_model)
  
if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_treefix datadir subst_model cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  cores = int(sys.argv[3])

  run_treefix_on_families(datadir, subst_model, cores)


