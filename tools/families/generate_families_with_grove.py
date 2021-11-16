import os
import sys
import shutil
import glob
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
from ete3 import SeqGroup
from ete3 import Tree


def get_presence_matrix_dict(presence_file):
  d = {}
  for line in open(presence_file).readlines()[1:]:
    line = line.replace("\n", "")
    sp = line.split(" ")
    d[sp[0]] = sp[1:]
  return d

def extract_from_grove_output(outputdir):
  generated_dirs = os.listdir(outputdir)
  assert(len(generated_dirs) == 1)
  grove_dir = os.path.join(outputdir, generated_dirs[0])
  fam.init_top_directories(outputdir)
  presence_matrix_dict = get_presence_matrix_dict(os.path.join(grove_dir, "iqt.pr_ab_matrix"))
  
  species_tree = os.path.join(grove_dir, "tree_best.newick")
  output_species_tree = fam.get_species_tree(outputdir)
  tree = Tree(species_tree)
  # arbitrarily root the tree  
  tree.set_outgroup(tree.get_leaves()[0])
  shutil.copy(species_tree, output_species_tree)
  genes_to_species = {}
  alignments = glob.glob(grove_dir + "/*.Fasta*")
  for taxon in presence_matrix_dict:
    assert(len(alignments) == len(presence_matrix_dict[taxon]))
  
  for ali in alignments:
    family = os.path.basename(ali).replace(".Fasta", "").replace(".", "_")
    print(family)
    family_index = int(family.split("part")[1])
    fam.init_family_directories(outputdir, family)
    output_ali = fam.get_alignment(outputdir, family)
    full_msa = SeqGroup(ali)
    pruned_msa = SeqGroup()
    species_to_genes = {}
    for entry in full_msa.get_entries():
      gene = entry[0]
      seq = entry[1]
      if (presence_matrix_dict[gene][family_index] == "1"):
        species_to_genes[gene] = [gene]
        pruned_msa.set_seq(gene, seq)
    pruned_msa.write("fasta", output_ali)
    mapping_file = fam.get_mappings(outputdir, family)
    fam.write_phyldog_mapping(species_to_genes, mapping_file)
  fam.postprocess_datadir(outputdir)
   
def run_grove(output_dir, minmiss, maxmiss, maxnodes, minfam, seed):
  command = []
  command.append(exp.raxmlgrove_exec)
  command.append("generate")
  command.append("-q")
  query = "MISSING_DATA_RATE > " + str(minmiss) + " and MISSING_DATA_RATE <= " + str(maxmiss) + " and MISSING_DATA_RATE != 'None' and NUM_TAXA < " + str(maxnodes) + " and RATE_AC != 'None' " + "and OVERALL_NUM_PARTITIONS >= " + str(minfam)
  command.append(query)
  command.append("--num-sequences")
  command.append("1")
  command.append("--insert-matrix-gaps")
  command.append("-o")
  command.append(output_dir)
  command.append("--seed")
  command.append(str(seed))
  print (" ".join(command))
  subprocess.check_call(command)

def get_output_dir(minmiss, maxmiss, maxnodes, minfam, seed):
  name = "grove"
  name += "_minmiss" + str(minmiss) 
  name += "_maxmiss" + str(maxmiss) 
  name += "_maxnodes" + str(maxnodes) 
  name += "_minfam" + str(minfam) 
  name += "_seed" + str(seed)
  output_dir = os.path.join(exp.families_datasets_root, name)
  return output_dir

def generate(minmiss, maxmiss, maxnodes, minfam, seed):
  output_dir = get_output_dir(minmiss, maxmiss, maxnodes, minfam, seed)
  os.mkdir(output_dir)
  exp.reset_dir(output_dir)
  run_grove(output_dir, minmiss, maxmiss, maxnodes, minfam, seed)
  extract_from_grove_output(output_dir)

if (__name__ == "__main__"): 
  if (len(sys.argv) < 5):
    print("Syntax: python" + os.path.basename(__file__) + " minmiss maxnodes minfam seed")
    sys.exit(1)
  minmiss = float(sys.argv[1])
  maxmiss = float(sys.argv[2])
  maxnodes = int(sys.argv[3])
  minfam = int(sys.argv[4])
  seed = int(sys.argv[5])
  generate(minmiss, maxmiss, maxnodes, minfam, seed)

