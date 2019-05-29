import os
import sys
import time
import launch_generax
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import experiments as exp
import saved_metrics

def run(datadir, starting_tree, with_transfer, cores, resultsdir):
  additional_arguments = []
  additional_arguments.append("--rec-model")
  dataset = os.path.basename(os.path.normpath(datadir))
  method = "GeneRax-"
  metric_name = "scaling-" + starting_tree + "-"
  if (with_transfer):
    metric_name += "DTL"
    additional_arguments.append("UndatedDTL")
  else:
    metric_name += "DL"
    additional_arguments.append("UndatedDL")
  start_time = time.time()
  launch_generax.run(dataset, "SPR", starting_tree, cores, additional_arguments, resultsdir, False)
  elapsed_time = time.time() - start_time
  saved_metrics.save_metrics(datadir, str(cores), str(elapsed_time), metric_name)

def launch(datadir, starting_tree, with_transfer, cluster, cores):
  command = ["python"]
  command.append(os.path.realpath(__file__))
  command.append(datadir)
  command.append(starting_tree)
  command.append(str(with_transfer))
  command.append(cluster)
  command.append(str(cores))
  command.append("--exprun")
  dataset = os.path.basename(os.path.normpath(datadir))
  print(command)
  resultsdir = os.path.join("ScalingGeneRax", dataset, starting_tree, "run")
  resultsdir = exp.create_result_dir(resultsdir, [])
  command.append(resultsdir)
  submit_path = os.path.join(resultsdir, "submit.sh")
  exp.submit(submit_path, " ".join(command), cores, cluster) 


if (__name__ == "__main__"):
  print(sys.argv)
  results_dir = exp.getAndDelete("--exprun", sys.argv, "")  
  print(sys.argv)
  is_run = (len(results_dir) > 0)
  if (len(sys.argv) != 6):
    print("Syntax: datadir starting_tree with_transfer cluster cores")
    exit(1)
  datadir = sys.argv[1]
  starting_tree = sys.argv[2]
  with_transfer = (int(sys.argv[3]) != 0)
  cluster = sys.argv[4]
  cores = int(sys.argv[5])
  if (is_run):
    run(datadir, starting_tree, with_transfer, cores, results_dir)
  else:
    launch(datadir, starting_tree, with_transfer, cluster, cores)

  
