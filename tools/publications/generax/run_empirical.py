import subprocess
import os
import sys
import common

datasets = []

datasets.append("cyano_simulated")

#common.run_all_reference_methods(datasets)
common.run_all_generax(datasets)

