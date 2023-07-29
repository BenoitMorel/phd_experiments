import sys
import os
import fam
from translate_gene_tree import translate  

def export(datadir, method, model, output):
  os.mkdir(output)
  for family in fam.get_families_list(datadir):
    gene_tree = fam.get_gene_tree_path(datadir, family, method, model)
    translation = translate(gene_tree, None)
    with open(os.path.join(output, family + "_" + method + "_" + model + ".newick"), "w") as writer:
        writer.write(translation)

if (__name__ == "__main__"): 
  if (len(sys.argv) != 5): 
    print("Syntax: python " + os.path.basename(__file__) + " datadir method subst_model outputdir")
    exit(1)

  datadir = sys.argv[1]
  method = sys.argv[2]
  model = sys.argv[3]
  output = sys.argv[4]
  export(datadir, method, model, output)


