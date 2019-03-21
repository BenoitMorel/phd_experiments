import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import subprocess

sys.path.insert(0, 'scripts')
import experiments as exp


def parameters_to_dataset(species, sites, bl, dup_rate):
  loss_rate = str(float(rates) / 2.0)
  return "jsim_s" + species + "_f50_sites" + sites + "_dna4_bl" + bl + "_d" + dup_rate + "_" + loss_rate

def get_param_from_dataset_name(parameter, dataset):
  if (parameter == "species"):
    return dataset.split("_")[1][1:]
  elif (parameter == "sites"):
    return dataset.split("_")[3][5:]
  elif (parameter == "bl"):
    return dataset.split("_")[5][2:]
  elif (parameter == "dup_rates"):
    return dataset.split("_")[6][1:]
  elif (parameter == "loss_rates"):
    return dataset.split("_")[7][1:]
  elif (parameter == "dl_ratio"):
    d = get_param_from_dataset_name("dup_rates", dataset)
    l = get_param_from_dataset_name("loss_rates", dataset)
    res =  str(float(d)/float(l))
    print(dataset + " " + res)
    return res
  else:
    return "invalid"

def get_available_datasets(prefix):
  res = []
  for dataset in os.listdir(exp.families_datasets_root):
    try:
      if (not dataset.startswith(prefix)):
        continue
      res.append(dataset)
    except:
      continue
  return res

def get_rf_from_logs(logs):
  lines = logs.split("\n")
  dico = {}
  for line in lines:
    if (line.startswith("- True - ")):
      line = line.replace("\t", " ")
      split = line.split(" ")
      method = split[3][:-1]
      rf = split[-1]
      dico[method] = rf
  return dico

def get_results(dataset):
  try:
    analyse_script = os.path.join(exp.tools_root, "families", "analyze_dataset.py")
    families_path = os.path.join(exp.benoit_datasets_root, "families", dataset, "families")
    cmd = []
    cmd.append("python")
    cmd.append(analyse_script)
    cmd.append(families_path)
    logs = subprocess.check_output(cmd)
    return get_rf_from_logs(logs) 
  except:
    print("Failed to get RF distances from dataset " + dataset)
    return None



def get_datasets_to_plot(datasets_rf_dico, fixed_params_dico):

  datasets_to_plot = []
  for dataset in datasets_rf_dico:
    ok = True
    for fixed_param in fixed_params_dico:
      dataset_param_value = get_param_from_dataset_name(fixed_param, dataset)
      if (dataset_param_value != fixed_params_dico[fixed_param]):
        ok = False
        continue
    if (ok):
      datasets_to_plot.append(dataset)
  return datasets_to_plot


def plot(datasets_rf_dico, x_param, fixed_params_dico, x_label, output):
  datasets_to_plot = get_datasets_to_plot(datasets_rf_dico, fixed_params_dico)
  datasets_to_plot.sort(key = lambda t: float(get_param_from_dataset_name(x_param, t)))
  df = pd.DataFrame()
  f, ax = plt.subplots(1)
 
  
  methods = ["RAxML-NG", "Notung", "Phyldog", "Treerecs", "GeneRax-Raxml", "GeneRax-Random"]
  fake_df = {}
  fake_df[x_param] = []
  for method in methods:
    fake_df[method] = []
  for dataset in datasets_to_plot:
    fake_df[x_param].append(float(get_param_from_dataset_name(x_param, dataset)))
    rf_dico = datasets_rf_dico[dataset]
    for method in rf_dico:
      if (not method in methods):
        print("Unknown method " + method)
        continue
      fake_df[method].append(float(rf_dico[method]))
  for elem in fake_df:
    df[elem] = fake_df[elem]
  
  for method in methods:
    style = "solid"
    plt.plot(x_param, method, data=df, marker='x', linestyle = style, linewidth=2, label = method)
    ax.set_ylim(bottom=0)
  plt.xlabel(x_label[x_param])
  plt.ylabel('RF distance')
  plt.legend()
  plt.savefig(output)
  print("Saving result in " + output)
  plt.close()

datasets = get_available_datasets("jsim")
datasets_rf_dico = {}
for dataset in datasets:
  res = get_results(dataset)
  if (res != None):
    datasets_rf_dico[dataset] = get_results(dataset)


x_label = {}
x_label["species"] = "Number of taxa in the species tree"
x_label["dup_rates"] = "Duplication rate"
x_label["bl"] = "Gene tree branch length multiplier"
x_label["dl_ratio"] = "Ratio between duplication and loss rates"
x_label["sites"] = "Number of sites"


params_value_dico_sites = {}
params_value_dico_sites["species"] = "25"
params_value_dico_sites["dup_rates"] = "1.0"
params_value_dico_sites["bl"] = "1.0"
params_value_dico_sites["dl_ratio"] = "2.0"
plot(datasets_rf_dico, "sites", params_value_dico_sites, x_label,  "sites.png")

params_value_dico_sites = {}
params_value_dico_sites["species"] = "25"
params_value_dico_sites["sites"] = "500"
params_value_dico_sites["bl"] = "1.0"
params_value_dico_sites["dl_ratio"] = "2.0"
plot(datasets_rf_dico, "dup_rates", params_value_dico_sites, x_label, "rates.png")

params_value_dico_sites = {}
params_value_dico_sites["species"] = "25"
params_value_dico_sites["sites"] = "500"
params_value_dico_sites["dup_rates"] = "0.5"
params_value_dico_sites["dl_ratio"] = "2.0"
plot(datasets_rf_dico, "bl", params_value_dico_sites, x_label, "bl.png")

params_value_dico_sites = {}
params_value_dico_sites["sites"] = "500"
params_value_dico_sites["dup_rates"] = "0.5"
params_value_dico_sites["bl"] = "1.0"
params_value_dico_sites["dl_ratio"] = "2.0"
plot(datasets_rf_dico, "species", params_value_dico_sites, x_label, "species.png")

params_value_dico_sites = {}
params_value_dico_sites["species"] = "25"
params_value_dico_sites["sites"] = "500"
params_value_dico_sites["dup_rates"] = "0.5"
params_value_dico_sites["bl"] = "1.0"
plot(datasets_rf_dico, "dl_ratio", params_value_dico_sites, x_label, "dl_ratio.png")





