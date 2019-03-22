import subprocess
import os
import sys
import common

datasets = []

datasets.append("jsim_s5_f100_sites500_dna4_bl1.0_d0.5_l0.25_t0.0")

#common.run_all_reference_methods(datasets)
common.run_all_generax(datasets)

