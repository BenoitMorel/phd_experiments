import run_raxml_supportvalues as raxml
import sys
import os

if __name__ == "__main__":

  if (len(sys.argv) < 3):
    print("syntax: python " + os.path.basename(__file__) + " pargenesdir outputfile")
    sys.exit(1)
  pargenesdir = sys.argv[1]
  outputfile = sys.argv[2]
  raxml.gather_likelihoods(pargenesdir, outputfile)

