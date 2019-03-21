import os
import sys
import common




datasets = []
#datasets.append("jsimdtl_s25_f50_sites500_dna4_bl1.0_d0.25_l0.25_t0.25")
#datasets.append("jsimdtl_s25_f50_sites500_dna4_bl1.0_d0.25_l0.25_t0.5")
#datasets.append("jsimdtl_s25_f50_sites500_dna4_bl1.0_d0.25_l0.25_t0.1")
#datasets.append("jsimdtl_s25_f50_sites500_dna4_bl1.0_d0.25_l0.25_t1.0")

common.generate_all_datasets(datasets)
common.run_all_reference_methods(datasets)
common.run_all_generax(datasets)




