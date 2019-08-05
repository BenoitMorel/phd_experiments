import os
import sys
sys.path.insert(0, 'scripts/generax')
import scaling_generax



def launch_dataset(dataset, cores_set, with_transfers, subst_model):
  for cores in cores_set:
    scaling_generax.launch(dataset, subst_model, "raxml-ng", with_transfers, "haswell", cores)
#    scaling_generax.launch(dataset, subst_model, "random", with_transfers, "haswell", cores)


#cores_set = [4, 8, 16, 32, 64, 128, 256, 512]
cores_set = [1024] #4, 8, 16, 32, 64, 128, 256, 512]
subst_model = "LG+G"
launch_dataset("../BenoitDatasets/families/cyano_empirical", cores_set, 1, subst_model)
#launch_dataset("../BenoitDatasets/families/cyano_empirical", cores_set, 0, subst_model)

