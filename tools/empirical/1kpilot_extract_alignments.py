"""
Rename 1kpilot gene alignements as family.msa
"""

import os
import sys
import shutil

def convert(raw_data_dir, output_ali_dir):
  os.mkdir(output_ali_dir)
  for f in os.listdir(raw_data_dir):
    family = f.split(".")[1]
    src = os.path.join(raw_data_dir, f)
    dest = os.path.join(output_ali_dir, family + ".msa")
    shutil.copyfile(src, dest)


if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python " + os.path.basename(__file__) + " raw_data_dir output_ali_dir")
    sys.exit(1)
  convert(sys.argv[1], sys.argv[2])


