import sys
sys.path.insert(0, 'scripts/generax')
sys.path.insert(0, 'tools/families')

import fam
import scaling_generax

def plot_scaling():
  datasets = ["cyano_empirical"]
  for dataset in datasets:
    datadir = fam.get_datadir(dataset)
    scaling_generax.plot_scaling_metric(datadir)

