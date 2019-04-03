import os
import sys
import common




datasets = []
def add_fixed_point():
  datasets.append("jsimdtl_s16_f100_sites500_dna4_bl1.0_d0.15_l0.15_t0.15")

def add_sites():
  datasets.append("jsimdtl_s16_f100_sites100_dna4_bl1.0_d0.15_l0.15_t0.15")
  datasets.append("jsimdtl_s16_f100_sites250_dna4_bl1.0_d0.15_l0.15_t0.15")
  datasets.append("jsimdtl_s16_f100_sites750_dna4_bl1.0_d0.15_l0.15_t0.15")
  datasets.append("jsimdtl_s16_f100_sites1000_dna4_bl1.0_d0.15_l0.15_t0.15")

def add_dtl_rates_scaler():
  datasets.append("jsimdtl_s16_f100_sites500_dna4_bl1.0_d0.02_l0.02_t0.02")
  datasets.append("jsimdtl_s16_f100_sites500_dna4_bl1.0_d0.08_l0.08_t0.08")
  datasets.append("jsimdtl_s16_f100_sites500_dna4_bl1.0_d0.25_l0.25_t0.25")

def add_transfers():
  datasets.append("jsimdtl_s16_f100_sites500_dna4_bl1.0_d0.15_l0.15_t0.05")
  datasets.append("jsimdtl_s16_f100_sites500_dna4_bl1.0_d0.15_l0.15_t0.25")
  datasets.append("jsimdtl_s16_f100_sites500_dna4_bl1.0_d0.15_l0.15_t0.4")

def add_species():
  datasets.append("jsimdtl_s5_f100_sites500_dna4_bl1.0_d0.15_l0.15_t0.15")
  datasets.append("jsimdtl_s10_f100_sites500_dna4_bl1.0_d0.15_l0.15_t0.15")
  datasets.append("jsimdtl_s12_f100_sites500_dna4_bl1.0_d0.15_l0.15_t0.15")
  datasets.append("jsimdtl_s19_f100_sites500_dna4_bl1.0_d0.15_l0.15_t0.15")
  datasets.append("jsimdtl_s27_f100_sites500_dna4_bl1.0_d0.15_l0.15_t0.15")
  datasets.append("jsimdtl_s41_f100_sites500_dna4_bl1.0_d0.15_l0.15_t0.15")

def add_bl():
  datasets.append("jsimdtl_s16_f100_sites500_dna4_bl0.01_d0.15_l0.15_t0.15")
  datasets.append("jsimdtl_s16_f100_sites500_dna4_bl0.05_d0.15_l0.15_t0.15")
  datasets.append("jsimdtl_s16_f100_sites500_dna4_bl0.1_d0.15_l0.15_t0.15")
  datasets.append("jsimdtl_s16_f100_sites500_dna4_bl0.2_d0.15_l0.15_t0.15")
  datasets.append("jsimdtl_s16_f100_sites500_dna4_bl0.5_d0.15_l0.15_t0.15")

def add_small():
  datasets.append("jsimdtl_s5_f5_sites50_dna4_bl1.0_d0.15_l0.15_t0.15")
  datasets.append("jsimdtl_s5_f10_sites50_dna4_bl1.0_d0.15_l0.15_t0.15")
  datasets.append("jsimdtl_s5_f15_sites50_dna4_bl1.0_d0.15_l0.15_t0.15")


add_fixed_point()
add_sites()
add_dtl_rates_scaler()
add_transfers()
add_bl()
add_species()
#add_small()

#common.generate_all_datasets(datasets)
#common.run_all_reference_methods(datasets)
#common.run_all_ALE(datasets, 1)
#common.run_all_generax(datasets)
common.run_all_analyzes(datasets)




