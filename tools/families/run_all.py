import os
import sys
import run_raxml_supportvalues as raxml
import run_treerecs as treerecs
import run_phyldog_light as phyldog


if __name__ == "__main__":
  if (len(sys.argv) != 4):
    print("syntax: python run_raxml_all.py dataset_dir bs_trees cores")
    sys.exit(1)
  dataset_dir = sys.argv[1]
  bs_trees = sys.argv[2]
  cores = int(sys.argv[3])
  raxml.run_pargenes_and_extract_trees(dataset_dir, bs_trees, cores)
  treerecs.run_treerecs_on_families(dataset_dir, cores)
  phyldog.run_phyldog_light_on_families(dataset_dir, cores)
  
  

