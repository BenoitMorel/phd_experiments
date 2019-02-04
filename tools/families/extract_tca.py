import sys
import os
sys.path.insert(0, 'tools/raxml/')
import raxml_get_tca_score as tca



def extract_tca_values(dataset_dir):
  concatenated_dir = os.path.join(dataset_dir, "pargenes", "concatenated_bootstraps")
  families_dir = os.path.join(dataset_dir, "families")
  for family in os.listdir(families_dir):
    tca_score = 0.0
    #try:
    bs_trees_file = os.path.join(concatenated_dir, family + ".bs")
    tca_score = tca.get_tca(bs_trees_file, family)
    #except:
      #print("Failed to get TCA score for " + family)
      #continue
    output = os.path.join(families_dir, family, "tca.txt")
    open(output, "w").writer.write(str(tca))

if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax: python script.py dataset_dir")
    exit(1)
  dataset_dir = sys.argv[1]
  extract_tca_values(dataset_dir)
