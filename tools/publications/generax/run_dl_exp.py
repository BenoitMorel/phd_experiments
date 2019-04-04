import subprocess
import os
import sys
import re
import common

datasets = []

def get_param_position(fixed_point, param_name):
  split = fixed_point.split("_")
  for i in range(0, len(split)):
    if (param_name ==  re.sub("[0-9]*[\.]*[0-9]*", "", split[i])):
      return i
  print("ERROR: unknown parameter " + param_name)
  exit(1)

def add_dataset(fixed_point, strings_to_replace):
  for string_to_replace in strings_to_replace:
    elems_to_replace = string_to_replace.split("_")
    split = fixed_point.split("_")
    for elem in elems_to_replace:
      param_name =  re.sub("[0-9]*[\.]*[0-9]*", "", elem)
      param_value =  re.sub("[a-zA-Z]*", "", elem)
      param_position = get_param_position(fixed_point, param_name)
      split[param_position] = param_name + str(param_value)
    dataset = "_".join(split)
    print("Add " + dataset)
    if (dataset in datasets):
      print("duplicate: " + dataset)
      exit(1)
    datasets.append(dataset)

if (False):
  datasets.append("jsim_s5_f10_sites200_dna4_bl0.5_d0.25_l0.25_t0.0")
  common.generate_all_datasets(datasets)
  common.run_all_reference_methods(datasets)
  common.run_all_generax(datasets)
  exit(0)

if (False):
  fixed_point_dl = "jsim_s19_f100_sites500_dna4_bl0.5_d0.25_l0.25_t0.0"
  datasets.append(fixed_point_dl)
  add_dataset(fixed_point_dl, ["s5", "s10", "s27", "s41"])
  add_dataset(fixed_point_dl, ["sites100", "sites250", "sites750", "sites1000"])
  add_dataset(fixed_point_dl, ["bl0.01", "bl0.05", "bl0.1", "bl0.2", "bl1.0", "bl2.0"])
  add_dataset(fixed_point_dl, ["f10", "f25", "f50", "f200"])
  add_dataset(fixed_point_dl, ["d0.01_l0.01", "d0.05_l0.05", "d0.1_l0.1", "d0.4_l0.4"])
  add_dataset(fixed_point_dl, ["d0.1", "d0.2", "d0.3", "d0.4"])
  common.generate_all_datasets(datasets)
  common.run_all_reference_methods(datasets)
  common.run_all_generax(datasets)

datasets = []

if (False):
  fixed_point_dtl = "jsimdtl_s16_f100_sites500_dna4_bl0.5_d0.1_l0.2_t0.1"
  datasets.append(fixed_point_dtl)
  add_dataset(fixed_point_dtl, ["s5", "s10", "s12", "s19", "s27", "s41"])
  add_dataset(fixed_point_dtl, ["sites100", "sites250", "sites750", "sites1000"])
  add_dataset(fixed_point_dtl, ["bl0.01", "bl0.05", "bl0.1", "bl0.2", "bl1.0", "bl2.0"])
  add_dataset(fixed_point_dtl, ["f10", "f25", "f50", "f200"])
  add_dataset(fixed_point_dtl, ["d0.01_l0.02_t0.01", "d0.05_l0.1_t0.05", "d0.15_l0.3_t0.15", "d0.2_l0.4_t0.2"])
  add_dataset(fixed_point_dtl, ["d0.01_l0.2_t0.19", "d0.05_l0.2_t0.15", "d0.15_l0.2_t0.05", "d0.19_l0.2_t0.01"])

  common.generate_all_datasets(datasets)
  common.run_all_reference_methods(datasets)
  common.run_all_generax(datasets)

datasets = []

if (False):
  datasets.append("cyano_simulated")
  #common.run_all_reference_methods(datasets)
  common.run_all_ALE(datasets, 0)
  common.run_all_generax(datasets)



#common.run_all_analyzes(datasets)

