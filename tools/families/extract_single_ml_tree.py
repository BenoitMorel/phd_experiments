import os
import sys
import fam

def extract(datadir, subst_model, index):
  for family in fam.get_families_list(datadir):
    f = fam.get_raxml_multiple_trees(datadir, subst_model, family)
    tree = open(f).readlines()[index]
    output = fam.build_gene_tree_path(datadir, subst_model, family, "raxml-ng-" + str(index))
    open(output, "w").write(tree)


if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python " + os.path.basename(__file__) + " datadir subst_model index")
    sys.exit(1)
  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  index = int(sys.argv[3])
  extract(datadir, subst_model, index)

