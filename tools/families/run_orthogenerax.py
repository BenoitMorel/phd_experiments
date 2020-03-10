import sys
import os
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'scripts/generax')
sys.path.insert(0, 'tools/family')
import launch_generax
import experiments as exp
import run_concatenation
import fam
import shutil
import time
import saved_metrics

def run_generax(datadir, subst_model, species_method, cores):
  dataset = os.path.basename(os.path.normpath(datadir))
  strategy = "SPR"
  species_tree = species_method
  starting_tree = "random"
  additional_arguments = []
  additional_arguments.append("--build-supermatrix")
  resultsdir = os.path.join("OrthoGeneRax", dataset, subst_model, species_tree)
  resultsdir = exp.create_result_dir(resultsdir, additional_arguments)
  return launch_generax.run(dataset, subst_model, strategy, species_tree, starting_tree, cores, additional_arguments, resultsdir)

def run_raxml_ng(datadir, subst_model, generax_output_dir, prefix, cores):
  supermatrix_path = os.path.join(generax_output_dir,"generax",  "superMatrix.fasta")
  partition_path = supermatrix_path + ".part"
  run_concatenation.run_raxml(subst_model, cores, prefix, supermatrix_path, partition_path)
  


def run_orthogenerax(datadir, subst_model, species_method, cores):
  method_name = "orthogenerax-" + species_method
  run_dir = fam.get_run_dir(datadir, subst_model, method_name)
  shutil.rmtree(run_dir, True)
  os.makedirs(run_dir)
  start = time.time()
  generax_output_dir = run_generax(datadir, subst_model, species_method, cores)
  run_raxml_ng(datadir, subst_model, generax_output_dir, run_dir, cores)
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name(method_name, subst_model), time1, "runtimes") 
  src = os.path.join(run_dir, "concatenation.raxml.bestTree")
  dest = fam.get_species_tree(datadir, subst_model, method_name)
  shutil.copy(src, dest)

if (__name__ == "__main__"): 
  if (len(sys.argv) != 5):
    print("Syntax: python " + os.path.basename(__file__) + " datadir subst_model species_method cores")
    sys.exit(1)
  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  species_method = sys.argv[3]
  cores = int(sys.argv[4])
  run_orthogenerax(datadir, subst_model, species_method, cores)

