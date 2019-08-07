import ete3
import sys

  
random_alignment_format = "fasta"

def create_random_tree(msa_file, output_file):
  global random_alignment_format
  msa = None
  msa = ete3.SeqGroup(msa_file, format="phylip_relaxed")

  tree = ete3.Tree()
  tree.populate(len(msa.get_entries()))

  leaves = tree.get_leaves()
  index = 0

  for entry in msa.get_entries():
    leaves[index].add_feature("name", entry[0])
    index += 1

  tree.write(outfile=output_file, format=1)
  

if (__name__ == "__main__"): 
  if (len(sys.argv) != 3):
    print("Syntax: python create_random_tree.py msa output_random_tree")
    sys.exit(1)

  msa_file = sys.argv[1]
  output_file = sys.argv[2]
  create_random_tree(msa_file, output_file)
