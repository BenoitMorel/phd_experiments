import os
import sys
import fam

if (__name__== "__main__"):
  if len(sys.argv) != 2:
    print("Syntax error: python " + os.path.basename(__file__) + " datadir")
    sys.exit(0)
  fam.postprocess_datadir(sys.argv[1])
  
