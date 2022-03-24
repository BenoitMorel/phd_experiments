import sys
import os
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/msa_edition')
import fam
sys.path.insert(0, 'scripts')
import experiments as exp
import glob
import re
import ete3
import read_msa

def get_genes(msa_path):
  msa = read_msa.read_msa(msa_path)
  genes = set()
  for entry in msa.get_entries():
    genes.add(entry[0])
  return genes

def copy_gene_trees_from(input_datadir, output_datadir, gene_tree_pattern, prune_tree):
  count = 0
  for family in fam.get_families_list(output_datadir):
    print(family)
    input_tree_dir = fam.get_gene_tree_dir(input_datadir, family)
    output_tree_dir = fam.get_gene_tree_dir(output_datadir, family)
    genes = {}
    if (prune_tree):
      genes = get_genes(fam.get_alignment(output_datadir, family))
    for gene_tree in os.listdir(input_tree_dir):
      if (gene_tree_pattern.match(gene_tree)):
        count += 1
        src = os.path.join(input_tree_dir, gene_tree)
        dest = os.path.join(output_tree_dir, gene_tree)
        if (prune_tree):
          tree = ete3.Tree(src)
          print(genes)
          tree.prune(genes)
          tree.write(outfile = dest, format = 1)
        else:
          exp.relative_symlink(src, dest)
  print("Copied " + str(count) + " gene trees")


if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax python " + os.path.basename(__file__) + " input_datadir output_datadir_pattern gene_tree_pattern prune_tree")
    sys.exit(1)

  input_datadir = sys.argv[1]
  output_datadir_pattern = sys.argv[2]
  gene_tree_pattern = re.compile(sys.argv[3])
  prune_tree = sys.argv[4] == "1"
  for output_datadir in glob.glob(output_datadir_pattern):
    print("Treating " + output_datadir)
    copy_gene_trees_from(input_datadir, output_datadir, gene_tree_pattern, prune_tree)




