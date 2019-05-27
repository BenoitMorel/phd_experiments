import os
import sys

import run_ALE


if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax: datadir cores")
    exit(1)
  datadir = sys.argv[1]
  cores = int(sys.argv[2])
  run_ALE.restart_exabayes_and_ALE(datadir, cores)
