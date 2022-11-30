
import os
import sys

sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "trees"))
import experiments as exp
import fam
import saved_metrics
import get_dico
import pairwise_rf_distance


def compare(datadir, gene_method, subst_model, cores):
  distances = pairwise_rf_distance.compute(datadir, gene_method, subst_model, cores)
  for family in fam.get_families_list(datadir):
    diff = fam.get_pythia_score(datadir, family)
    print(family + "\t" + str(diff) + "\t" + str(distances[family]))


if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax python " + os.path.basename(__file__) + " datadir gene_method subst_model 40")
    sys.exit(1)
  datadir = sys.argv[1]
  gene_method = sys.argv[2]
  subst_model = sys.argv[3]
  cores=  sys.argv[4]
  compare(datadir, gene_method, subst_model, cores)


