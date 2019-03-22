import subprocess
import os
import sys
import common

datasets = []
# fix point:
#datasets.append("jsim_s12_f50_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")


# bl
#datasets.append("jsim_s12_f50_sites500_dna4_bl0.25_d0.5_l0.25_t0.0")
#datasets.append("jsim_s12_f50_sites500_dna4_bl0.5_d0.5_l0.25_t0.0")
#datasets.append("jsim_s12_f50_sites500_dna4_bl2.0_d0.5_l0.25_t0.0")
#datasets.append("jsim_s12_f50_sites500_dna4_bl3.0_d0.5_l0.25_t0.0")

# sites
#datasets.append("jsim_s12_f50_sites100_dna4_bl1.0_d0.5_l0.25_t0.0")
#datasets.append("jsim_s12_f50_sites250_dna4_bl1.0_d0.5_l0.25_t0.0")
#datasets.append("jsim_s12_f50_sites750_dna4_bl1.0_d0.5_l0.25_t0.0")
#datasets.append("jsim_s12_f50_sites1000_dna4_bl1.0_d0.5_l0.25_t0.0")

# rates multiplier
#datasets.append("jsim_s12_f50_sites500_dna4_bl1.0_d0.1_l0.05_t0.0")
#datasets.append("jsim_s12_f50_sites500_dna4_bl1.0_d0.2_l0.1_t0.0")
#datasets.append("jsim_s12_f50_sites500_dna4_bl1.0_d1_l0.5_t0.0")
#datasets.append("jsim_s12_f50_sites500_dna4_bl1.0_d2_l1.0_t0.0")

# dl ratio
#datasets.append("jsim_s12_f50_sites500_dna4_bl1.0_d0.5_l0.05_t0.0")
#datasets.append("jsim_s12_f50_sites500_dna4_bl1.0_d0.5_l0.1_t0.0")
#datasets.append("jsim_s12_f50_sites500_dna4_bl1.0_d0.5_l0.5_t0.0")
#datasets.append("jsim_s12_f50_sites500_dna4_bl1.0_d0.5_l1.0_t0.0")

# species:
#datasets.append("jsim_s8_f50_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")
#datasets.append("jsim_s13_f50_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")
#datasets.append("jsim_s16_f50_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")
#datasets.append("jsim_s21_f50_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")
#datasets.append("jsim_s34_f50_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")

# species ++:
datasets.append("jsim_s5_f100_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")
datasets.append("jsim_s10_f100_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")
datasets.append("jsim_s19_f100_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")
datasets.append("jsim_s27_f100_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")
datasets.append("jsim_s41_f100_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")


common.generate_all_datasets(datasets)
common.run_all_reference_methods(datasets)
common.run_all_generax(datasets)

