import ete3
import sys

if (len(sys.argv) != 3):
  print("Syntax: python create_random_tree.py msa output_random_tree")
  sys.exit(1)

msa_file = sys.argv[1]
output_file = sys.argv[2]

msa = None
try:
  msa = ete3.SeqGroup(msa_file, format="fasta")
except:
  msa = ete3.SeqGroup(msa_file, format="phylip")

tree = ete3.Tree()
tree.populate(len(msa.get_entries()))

leaves = tree.get_leaves()
index = 0

for entry in msa.get_entries():
  leaves[index].add_feature("name", entry[0])
  index += 1

tree.write(outfile=output_file, format=1)

