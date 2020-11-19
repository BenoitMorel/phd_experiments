import os
import sys
import fam_data


# get a list of replicated datasets with varying parameters
# from a reference dataset, a list of parameters to vary, 
# and a range of replicates
def get_dataset_list(ref_dataset, strings_to_replace, replicates):
  seeds_to_replace = []
  for i in replicates:
    seeds_to_replace.append("seed" + str(i))
  unreplicated_datasets = []
  fam_data.get_dataset_variations(unreplicated_datasets, ref_dataset, strings_to_replace)
  replicated_datasets = []
  for d in unreplicated_datasets:
    fam_data.get_dataset_variations(replicated_datasets, d, seeds_to_replace)
  return replicated_datasets

