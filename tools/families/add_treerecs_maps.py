import sys
import os
import shutil
sys.path.insert(0, 'tools/treerecs')
import phyldog_to_treerecs_map


if (len(sys.argv) != 2):
  print("Syntax: python script.py dataset_dir")

dataset_dir = sys.argv[1]
families_dir = os.path.join(dataset_dir, "families")
for family in os.listdir(families_dir):
  phyldog_mapping = os.path.join(families_dir, family, "mapping.link")
  treerecs_mapping = os.path.join(families_dir, family, "treerecs_mapping.link")
  try:
    phyldog_to_treerecs_map.convert(phyldog_mapping, treerecs_mapping)
  except:
    print("error")

