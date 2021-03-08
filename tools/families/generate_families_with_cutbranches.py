import sys
import os
import fam
sys.path.insert(0, 'tools/trees')
import cut_long_branches
import ete3



ali = "../BenoitDatasets/families/pdb_plants23/families/Phy003MBZY_CUCME/alignment.msa"
tree= "../BenoitDatasets/families/pdb_plants23/families/Phy003MBZY_CUCME/gene_trees/raxml-ng.bestAA.geneTree.newick"

out = "plop"


subtrees = cut_long_branches.cut_long_branches(tree, 2.0)
idx = 1
for subtree in subtrees:
  leaves = set(subtree.get_leaf_names())
  if (len(leaves) < 4):
    continue
  newtree = os.path.join(out, "tree." + str(idx) + ".newick")
  newali = os.path.join(out, "ali." + str(idx) + ".fasta")
  idx = idx + 1
  subtree.write(outfile = newtree)
  seqs = ete3.SeqGroup(ali, format="fasta")
  newseqs = ete3.SeqGroup()
  for seq in  seqs.get_entries():
    if (seq[0] in leaves):
      newseqs.set_seq(seq[0], seq[1])

  newseqs.write(outfile = newali, format = "fasta")

