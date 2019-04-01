import subprocess
import os
import sys
import common

datasets = []

datasets.append("cyano_simulated")
#datasets.append("bench_libpll")

#common.run_all_reference_methods(datasets)
common.run_all_raxml_light(datasets)
common.run_all_generax(datasets)

