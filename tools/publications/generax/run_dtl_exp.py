import os
import sys
import common




datasets = []
# fixed point
#datasets.append("jsimdtl_s12_f50_sites500_dna4_bl1.0_d0.25_l0.25_t0.25")


# sites
#datasets.append("jsimdtl_s12_f50_sites100_dna4_bl1.0_d0.25_l0.25_t0.25")
#datasets.append("jsimdtl_s12_f50_sites250_dna4_bl1.0_d0.25_l0.25_t0.25")
#datasets.append("jsimdtl_s12_f50_sites750_dna4_bl1.0_d0.25_l0.25_t0.25")
#datasets.append("jsimdtl_s12_f50_sites1000_dna4_bl1.0_d0.25_l0.25_t0.25")

# DTL rates scaler
datasets.append("jsimdtl_s12_f50_sites500_dna4_bl1.0_d0.01_l0.01_t0.01")
datasets.append("jsimdtl_s12_f50_sites500_dna4_bl1.0_d0.1_l0.1_t0.1")
datasets.append("jsimdtl_s12_f50_sites500_dna4_bl1.0_d0.5_l0.5_t0.5")
datasets.append("jsimdtl_s12_f50_sites500_dna4_bl1.0_d1.0_l1.0_t1.0")


#transfers
#datasets.append("jsimdtl_s12_f50_sites500_dna4_bl1.0_d0.25_l0.25_t0.5")
#datasets.append("jsimdtl_s12_f50_sites500_dna4_bl1.0_d0.25_l0.25_t0.1")
#datasets.append("jsimdtl_s12_f50_sites500_dna4_bl1.0_d0.25_l0.25_t1.0")

common.generate_all_datasets(datasets)
common.run_all_reference_methods(datasets)
common.run_all_generax(datasets)




