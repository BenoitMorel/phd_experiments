import os
import sys
sys.path.insert(0, 'scripts/generax')
import scaling_generax



def launch_dataset(dataset, cores_set):
  for cores in cores_set:
    #scaling_generax.launch(dataset, "raxml-ng", 0, "haswell", cores)
    scaling_generax.launch(dataset, "raxml-ng", 1, "haswell", cores)


cores_set = [8, 16, 32, 64, 128, 256, 512, 1024]
launch_dataset("../BenoitDatasets/families/ensembl_96_ncrna_primates", cores_set) 
