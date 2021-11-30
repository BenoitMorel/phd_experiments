import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "families"))
import experiments as exp
import create_random_tree

def simulate_msa(msa_path, tree_path, sites, outputdir): 
  
  #gene_tree = os.path.join(genes_dir, "gene_" + str(i) + ".newick")
  
  #seqgen_gene_tree = os.path.join(genseq_genes_dir, "gene_" + str(i) + ".newick")
  #shutil.copyfile(gene_tree, seqgen_gene_tree)
  #subprocess.check_call(["sed", "-i", "s/[SDLT]//g", seqgen_gene_tree])
  command = []
  command.append(exp.seq_gen_exec)
  command.append("-l")
  command.append(str(sites))
  command.append("-m")
  command.append("GTR")
  command.append("-of")
  command.append(tree_path)
  #command.append("-z")
  #command.append(str(int(i) + int(seed)))
  with open(msa_path, "w") as writer:
    subprocess.check_call(command, stdout=writer)
  print("Simulated MSA saved in " + msa_path)

def simulate(taxa, sites, outputdir):
  os.mkdir(outputdir)  
  tree = create_random_tree.create_random_tree_taxa_number(taxa)
  tree_path = os.path.join(outputdir, "simulated_tree.newick")
  msa_path = os.path.join(outputdir, "simulated_msa.fasta")
  tree.write(outfile = tree_path, format = 1)
  print("Simulated tree saved in " + tree_path)
  simulate_msa(msa_path, tree_path, sites, outputdir) 


if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python " + os.path.basename(__file__) + "taxa sites output_dir")
    sys.exit(1)
  taxa = int(sys.argv[1])
  sites = int(sys.argv[2])
  outputdir = sys.argv[3]
  simulate(taxa, sites, outputdir)
