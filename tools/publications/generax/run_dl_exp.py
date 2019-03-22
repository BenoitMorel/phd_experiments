import subprocess
import os
import sys
import common

datasets = []


def add_fixed_point():
  datasets.append("jsim_s19_f100_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")


def add_bl():
  datasets.append("jsim_s19_f100_sites500_dna4_bl0.25_d0.5_l0.25_t0.0")
  datasets.append("jsim_s19_f100_sites500_dna4_bl0.5_d0.5_l0.25_t0.0")
  datasets.append("jsim_s19_f100_sites500_dna4_bl2.0_d0.5_l0.25_t0.0")
  datasets.append("jsim_s19_f100_sites500_dna4_bl3.0_d0.5_l0.25_t0.0")

def add_sites():
  datasets.append("jsim_s19_f100_sites100_dna4_bl1.0_d0.5_l0.25_t0.0")
  datasets.append("jsim_s19_f100_sites250_dna4_bl1.0_d0.5_l0.25_t0.0")
  datasets.append("jsim_s19_f100_sites750_dna4_bl1.0_d0.5_l0.25_t0.0")
  datasets.append("jsim_s19_f100_sites1000_dna4_bl1.0_d0.5_l0.25_t0.0")

def add_rates_multiplier():
  datasets.append("jsim_s19_f100_sites500_dna4_bl1.0_d0.1_l0.05_t0.0")
  datasets.append("jsim_s19_f100_sites500_dna4_bl1.0_d0.2_l0.1_t0.0")
  datasets.append("jsim_s19_f100_sites500_dna4_bl1.0_d1_l0.5_t0.0")
  datasets.append("jsim_s19_f100_sites500_dna4_bl1.0_d2_l1.0_t0.0")

def add_dl_ratio():
  datasets.append("jsim_s19_f100_sites500_dna4_bl1.0_d0.5_l0.05_t0.0")
  datasets.append("jsim_s19_f100_sites500_dna4_bl1.0_d0.5_l0.1_t0.0")
  datasets.append("jsim_s19_f100_sites500_dna4_bl1.0_d0.5_l0.5_t0.0")
  datasets.append("jsim_s19_f100_sites500_dna4_bl1.0_d0.5_l1.0_t0.0")

def add_species():
  datasets.append("jsim_s5_f100_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")
  datasets.append("jsim_s10_f100_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")
  datasets.append("jsim_s27_f100_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")
  datasets.append("jsim_s41_f100_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")


#add_fixed_point()
add_bl()
add_sites()
add_rates_multiplier()
add_dl_ratio()
#add_species()

#common.generate_all_datasets(datasets)
common.run_all_reference_methods(datasets)
common.run_all_generax(datasets)

