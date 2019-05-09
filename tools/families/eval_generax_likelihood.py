import sys
import os
import shutil
sys.path.insert(0, 'scripts/generax')
import launch_generax
import saved_metrics

try:
      from StringIO import StringIO
except ImportError:
      from io import StringIO
import saved_metrics

def eval_likelihood(dataset_dir, starting_tree, with_transfers, is_protein, cores):
  temp_dir = os.path.join(dataset_dir, "temp", starting_tree)
  shutil.rmtree(temp_dir, ignore_errors  = True)
  os.makedirs(temp_dir)
  generax_families_file = os.path.join(temp_dir, "families")
  launch_generax.build_generax_families_file(dataset_dir, starting_tree, is_protein, generax_families_file)
  additional_arguments = []
  if (with_transfers):
    additional_arguments.append("--rec-model")
    additional_arguments.append("UndatedDTL")
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
  
def eval_and_save_likelihood(dataset_dir, starting_tree, with_transfers, is_protein, cores):
  stats = eval_likelihood(dataset_dir, starting_tree, with_transfers, is_protein, cores)
  for metric in stats:
    metric_name = metric
    if (with_transfers):
      metric_name += "_DTL"
    else:
      metric_name += "_DL"
    saved_metrics.save_metrics(dataset_dir, starting_tree, stats[metric], metric_name)
     
if (__name__ == "__main__"): 
  if (len(sys.argv) != 6): 
     print("Syntax: python " + os.path.basename(__file__) + " dataset_dir starting_tree with_transfers is_protein ")
     exit(1)
  dataset_dir = sys.argv[1]
  starting_tree = sys.argv[2]
  with_transfers = int(sys.argv[3])
  is_protein =  int(sys.argv[4])
  cores =  int(sys.argv[4])
  eval_and_save_likelihood(dataset_dir, starting_tree, with_transfers, is_protein, cores)

