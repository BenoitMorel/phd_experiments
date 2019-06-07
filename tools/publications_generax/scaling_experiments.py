import os
import sys
sys.path.insert(0, 'scripts/generax')
import scaling_generax



def launch_dataset(dataset, cores_set, with_transfers):
  for cores in cores_set:
    scaling_generax.launch(dataset, "raxml-ng", with_transfers, "haswell", cores)


cores_set = [8, 16, 32, 64, 128, 256, 512]
launch_dataset("../BenoitDatasets/families/ensembl_96_ncrna_primates", cores_set, 0) 
launch_dataset("../BenoitDatasets/families/cyano_empirical", cores_set, 1)

