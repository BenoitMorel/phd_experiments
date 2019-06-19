import sys
import os
import shutil
sys.path.insert(0, 'scripts/generax')
sys.path.insert(0, 'scripts')
import experiments as exp
import launch_generax
import saved_metrics
import fam

try:
      from StringIO import StringIO
except ImportError:
      from io import StringIO
import saved_metrics

def eval_likelihood(dataset_dir, run_name, with_transfers, subst_model, cores):
  temp_dir = os.path.join(dataset_dir, "temp", run_name)
  shutil.rmtree(temp_dir, ignore_errors  = True)
  os.makedirs(temp_dir)
  generax_families_file = os.path.join(temp_dir, "families")
  launch_generax.build_generax_families_file(dataset_dir, run_name, subst_model, generax_families_file)
  additional_arguments = []
  if (with_transfers):
    additional_arguments.append("--rec-model")
    additional_arguments.append("UndatedDTL")
  else:
    additional_arguments.append("--rec-model")
    additional_arguments.append("UndatedDL")
  sys.stdout = open(os.devnull, 'w')
  launch_generax.run_generax(dataset_dir, "EVAL", generax_families_file, "normal", cores, additional_arguments, temp_dir)
  sys.stdout = sys.__stdout__
  
  lines = open(os.path.join(temp_dir, "generax", "stats.txt")).readlines()
  stats = {}
  for line in lines:
    split = line.split()
    key = split[0].replace(":", "")
    value = split[1]
    stats[key] = value
  return stats
  
def eval_and_save_likelihood(dataset_dir, run_name, with_transfers, subst_model, cores):
  print("Evaluating likelihood for " + run_name + " and " + dataset_dir + " transfers=" + str(with_transfers))
  if (run_name == "all"):
    for run in fam.get_successful_runs(dataset_dir):
      eval_and_save_likelihood(dataset_dir, run, with_transfers, subst_model, cores)
    return

  stats = eval_likelihood(dataset_dir, run_name, with_transfers, subst_model, cores)
  for metric in stats:
    metric_name = metric
    if (with_transfers):
      metric_name += "_DTL"
    else:
      metric_name += "_DL"
    saved_metrics.save_metrics(dataset_dir, run_name, stats[metric], metric_name)

def launch(dataset_dir, run_name, subst_model, cluster, cores):
  dataset = os.path.basename(os.path.normpath(dataset_dir)) 
  command = ["python"]
  command.append(dataset_dir)
  command.append(run_name)
  command.append(str(subst_model))
  command.append(str(cores))
  command.append("--exprun")
  resultsdir = os.path.join("EvalGeneraxLikelihood", dataset, run_name)
  resultsdir = exp.create_result_dir(resultsdir, [])
  submit_path = os.path.join(resultsdir, "submit.sh")
  exp.submit(submit_path, " ".join(command), cores, cluster) 

     
if (__name__ == "__main__"): 
  is_run = ("--exprun" in sys.argv)
  if (len(sys.argv) < 6): 
     print("Syntax: python " + os.path.basename(__file__) + " dataset_dir run_name subst_model cores cluster")
     print("run_name can be all")
     exit(1)
  dataset_dir = sys.argv[1]
  run_name = sys.argv[2]
  subst_model =  int(sys.argv[3])
  cores =  int(sys.argv[4])
  if (is_run):
    eval_and_save_likelihood(dataset_dir, run_name, False, subst_model, cores)
    eval_and_save_likelihood(dataset_dir, run_name, True, subst_model, cores)
  else:
    cluster = sys.argv[5]
    launch(dataset_dir, run_name, subst_model, cluster, cores)
