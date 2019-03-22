import os
import sys
import common




datasets = []
def add_fixed_point():
  datasets.append("jsimdtl_s19_f100_sites500_dna4_bl1.0_d0.25_l0.25_t0.25")

def add_sites():
  datasets.append("jsimdtl_s19_f100_sites100_dna4_bl1.0_d0.25_l0.25_t0.25")
  datasets.append("jsimdtl_s19_f100_sites250_dna4_bl1.0_d0.25_l0.25_t0.25")
  datasets.append("jsimdtl_s19_f100_sites750_dna4_bl1.0_d0.25_l0.25_t0.25")
  datasets.append("jsimdtl_s19_f100_sites1000_dna4_bl1.0_d0.25_l0.25_t0.25")

def add_dtl_rates_scaler():
  datasets.append("jsimdtl_s19_f100_sites500_dna4_bl1.0_d0.01_l0.01_t0.01")
  datasets.append("jsimdtl_s19_f100_sites500_dna4_bl1.0_d0.1_l0.1_t0.1")
  datasets.append("jsimdtl_s19_f100_sites500_dna4_bl1.0_d0.5_l0.5_t0.5")
  datasets.append("jsimdtl_s19_f100_sites500_dna4_bl1.0_d1.0_l1.0_t1.0")

def add_transfers():
  datasets.append("jsimdtl_s19_f100_sites500_dna4_bl1.0_d0.25_l0.25_t0.5")
  datasets.append("jsimdtl_s19_f100_sites500_dna4_bl1.0_d0.25_l0.25_t0.1")
  datasets.append("jsimdtl_s19_f100_sites500_dna4_bl1.0_d0.25_l0.25_t1.0")

def add_species():
  datasets.append("jsimdtl_s5_f100_sites500_dna4_bl1.0_d0.25_l0.25_t0.25")
  datasets.append("jsimdtl_s10_f100_sites500_dna4_bl1.0_d0.25_l0.25_t0.25")
  datasets.append("jsimdtl_s27_f100_sites500_dna4_bl1.0_d0.25_l0.25_t0.25")
  datasets.append("jsimdtl_s41_f100_sites500_dna4_bl1.0_d0.25_l0.25_t0.25")

add_fixed_point()
add_sites()
add_dtl_rates_scaler()
add_transfers()

common.generate_all_datasets(datasets)
common.run_all_reference_methods(datasets)
common.run_all_generax(datasets)




