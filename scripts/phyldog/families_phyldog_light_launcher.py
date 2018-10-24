import os
import sys
sys.path.insert(0, 'scripts')
import experiments as exp



def generate_scheduler_commands_file(dataset_dir, cores, output_dir, scheduler_output_dir):
  families_dir = os.path.join(dataset_dir, "families")
  results_dir = os.path.join(scheduler_output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  

  with open(scheduler_commands_file, "w") as writer:
    for family in os.listdir(families_dir):
      family_dir = os.path.join(families_dir, family)
      species_tree = os.path.join(family_dir, "speciesTree.newick")
      phyldog_dir = os.path.join(family_dir, "phyldog")
      try:
        os.makedirs(phyldog_dir)
      except:
        pass
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append("species.tree.file=" + os.path.join(dataset_dir, "phyldogSpeciesTree.newick"))
      command.append("gene.tree.file=" + os.path.join(family_dir, "raxmlGeneTree.newick"))
      command.append("input.sequence.file=" + os.path.join(family_dir, "alignment.msa"))
      command.append("taxaseq.file=" + os.path.join(family_dir, "phyldogMapping.link"))
      command.append("likelihood.evaluator=LIBPLL2")
      command.append("model=GTR ")
      os.makedirs(os.path.join(results_dir, family))
      command.append("output.file=" + os.path.join(family_dir, "phyldog", "phyldog"))
      #command.append("output.file=" + os.path.join(results_dir, family, "phyldog"))
      #print(" ".join(command))
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def generate_scheduler_command(command_file, cores, scheduler_output_dir):
  command = ""
  parallelization = "onecore"
  command += "mpirun -np " + str(cores) + " "
  command += exp.mpischeduler_exec + " "
  command += "--" + parallelization + "-scheduler "
  command += exp.phyldog_light_exec + " "
  command += command_file + " "
  command += scheduler_output_dir + " " 
  command += "0 "
  return command 



max_args_number = 3
if len(sys.argv) < max_args_number:
  print("Syntax error: python families_phyldog_light_launcher.py cluster cores.")
  print("Cluster can be either normal, haswell or magny")
  sys.exit(0)


cluster = sys.argv[1]
cores = int(sys.argv[2])


resultsdir = os.path.join("phyldog_mpischeduler", cluster + "_" + str(cores), "run")
resultsdir = exp.create_result_dir(resultsdir)
result_msg = ""
exp.write_results_info(resultsdir, result_msg) 
output_dir = resultsdir 

datadir = os.path.join(exp.bigdatasets_root, "simuls_higher_rate")


scheduler_output_dir = os.path.join(output_dir, "scheduler_run")
os.makedirs(scheduler_output_dir)

scheduler_commands_file = generate_scheduler_commands_file(datadir, cores, output_dir, scheduler_output_dir)


command = generate_scheduler_command(scheduler_commands_file, cores,  scheduler_output_dir)


submit_path = os.path.join(resultsdir, "phyldog_scheduler_submit.sh")
print("Results will be in " + resultsdir)
exp.submit(submit_path, command, cores, cluster) 




