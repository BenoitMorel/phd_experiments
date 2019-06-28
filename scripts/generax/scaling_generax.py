import os
import sys
import time
import launch_generax
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/plotters')
import experiments as exp
import saved_metrics
import plot_line

def get_run_name(starting_tree, with_transfer, subst_model, cores):
  run_name = "scalinggenerax-" + starting_tree + "."
  if (with_transfer):
    run_name += "dtl."
  else:
    run_name += "dl."
  run_name += subst_model
  run_name += "-c" + str(cores)
  return run_name

def run(datadir, subst_model, starting_tree, with_transfer, cores, resultsdir):
  
  additional_arguments = []
  additional_arguments.append("--run")
  additional_arguments.append(get_run_name(starting_tree, with_transfer, subst_model, cores)) 
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
  launch_generax.run(dataset, subst_model, "SPR", starting_tree, cores, additional_arguments, resultsdir, False)
  elapsed_time = time.time() - start_time
  saved_metrics.save_metrics(datadir, str(cores), str(elapsed_time), metric_name)

def launch(datadir, subst_model, starting_tree, with_transfer, cluster, cores):
  command = ["python"]
  command.append(os.path.realpath(__file__))
  command.append(datadir)
  command.append(subst_model)
  command.append(starting_tree)
  command.append(str(with_transfer))
  command.append(cluster)
  command.append(str(cores))
  command.append("--exprun")
  dataset = os.path.basename(os.path.normpath(datadir))
  print(command)
  resultsdir = os.path.join("ScalingGeneRax", dataset, starting_tree, get_run_name(starting_tree, with_transfer, subst_model, cores))
  resultsdir = exp.create_result_dir(resultsdir, [])
  command.append(resultsdir)
  submit_path = os.path.join(resultsdir, "submit.sh")
  exp.submit(submit_path, " ".join(command), cores, cluster) 


def plot_scaling_metric(datadir):
  dataset = os.path.basename(os.path.normpath(datadir))
  print("dataset " + dataset)
  for metric_name in saved_metrics.get_all_metric_names(datadir):
    if (not metric_name.startswith("scaling-")):
      continue
    metrics = saved_metrics.get_metrics(datadir, metric_name)
    output_file = dataset + "-" + metric_name + ".svg"
    xvalues = []
    yvalues = []
    for key in metrics:
      xvalues.append(int(key))
      yvalues.append(float(metrics[key]))
    xvalues, yvalues = (list(t) for t in zip(*sorted(zip(xvalues, yvalues))))    
    print(xvalues)
    print(yvalues)
    cost_per_cores = yvalues[1] * xvalues[1]
    for i in range(0, len(xvalues)):
        yvalues[i] = cost_per_cores / yvalues[i]
    yvalues_list = [yvalues, xvalues]
    lines_captions = ["GeneRax", "Theoretical optimum"]
    plot_line.plot_line(xvalues, yvalues_list, "GeneRax parallel efficiency",  "Cores", "Speedup", output_file, lines_captions)
    print("Output file: " + output_file)

if (__name__ == "__main__"):
  results_dir = exp.getAndDelete("--exprun", sys.argv, "")  
  is_run = (len(results_dir) > 0)
  if (len(sys.argv) < 7):
    print("Syntax: datadir subst_model starting_tree with_transfer cluster cores")
    exit(1)
  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  starting_tree = sys.argv[3]
  with_transfer = (int(sys.argv[4]) != 0)
  cluster = sys.argv[5]
  cores = int(sys.argv[6])
  if (is_run):
    run(datadir, subst_model, starting_tree, with_transfer, cores, results_dir)
  else:
    launch(datadir, subst_model, starting_tree, with_transfer, cluster, cores)

  
