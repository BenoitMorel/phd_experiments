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
import ete3

def get_csv(output_dir, family):
  return os.path.join(output_dir, "train_" + family + ".csv")

def generate_scheduler_command_file(datadir, subst_model, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_dir = fam.get_family_path(datadir, family)
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      #../FastTree/FastTree -out out.newick -nt -gtr
      command.append(fam.get_alignment(datadir, family))
      command.append(model)
      command.append("1")
      command.append(get_csv(output_dir, family))
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def extract_csv(datadir, output_dir, model):
  output = os.path.join(fam.get_misc_dir(datadir), "train_" + model.replace("+", ""))
  output_ok = output + "_ok.csv"
  output_ll = output + "_ll.csv"
  write_ok = open(output_ok, "w")
  write_ll = open(output_ll, "w")
  for family in fam.get_families_list(datadir):
    for line in open(get_csv(output_dir, family), "r").readlines():
      line = line.replace("\n", "")
      sp = line.split(",")
      line_ok = sp[0:-2] + [sp[-1]]
      line_ll = sp[0:-2] + [sp[-2]]
      write_ok.write(",".join(line_ok) + "\n") 
      write_ll.write(",".join(line_ll) + "\n") 
      
  print("Trained data saved into:")
  print(output_ok)
  print(output_ll)

def train_deeprax(datadir, model, cores):          
  output_dir = fam.get_run_dir(datadir, model, "traindeeprax")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_command_file(datadir, model, cores, output_dir)
  
  exp.run_with_scheduler(exp.gensamples_exec, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  extract_csv(datadir, output_dir, model)
  

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python " + os.path.basename(__file__) + " datadir gene_trees subst_model")
    sys.exit(1)
  datadir = sys.argv[1]
  model = sys.argv[2]
  cores = int(sys.argv[3])
  train_deeprax(datadir, model, cores)          
  
  

