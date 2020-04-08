import sys
import os
sys.path.insert(0, 'tools/families')
import fam
sys.path.insert(0, 'scripts')
import experiments as exp
import glob

def copy_gene_trees_from(input_datadir, output_datadir):
  for family in fam.get_families_list(output_datadir):
    input_tree_dir = fam.get_gene_tree_dir(input_datadir, family)
    output_tree_dir = fam.get_gene_tree_dir(output_datadir, family)
    for gene_tree in os.listdir(input_tree_dir):
      src = os.path.join(input_tree_dir, gene_tree)
      dest = os.path.join(output_tree_dir, gene_tree)
      exp.relative_symlink(src, dest)



if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python " + os.path.basename(__file__) + " input_datadir output_datadir_pattern")
    sys.exit(1)

  input_datadir = sys.argv[1]
  output_datadir_pattern = sys.argv[2]
  for output_datadir in glob.glob(output_datadir_pattern):
    print("Treating " + output_datadir)
    copy_gene_trees_from(input_datadir, output_datadir)




