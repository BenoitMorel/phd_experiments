import os
import sys
import shutil

def find_family_index(alignments_info_file, alignment_name):
  lines = open(alignments_info_file).readlines()[1:]
  for i in range(0, len(lines)):
    if "/" + alignment_name in lines[i] or lines[i].startswith(alignment_name):
      return i
  return -1

def extract_line(input_file, index, output_file):
  lines = open(input_file).readlines()
  with open(output_file, "w") as writer:
    current_index = 0
    for line in lines:
      if (">" in line):
        continue
      else:
        if (current_index == index):
          writer.write(line)
          return
        else:
          current_index += 1

def extract_gene_family(input_dataset_dir, alignment_name, output):
  species_tree = os.path.join(input_dataset_dir, "speciesTree.newick")
  true_trees = os.path.join(input_dataset_dir, "trueGeneTrees.newick")
  raxml_trees = os.path.join(input_dataset_dir, "geneTrees.newick")
  treerecs_trees = os.path.join(input_dataset_dir, "treerecs_output.newick.best")
  alignments_info_file = os.path.join(input_dataset_dir, "alignment.txt")
  msa_file = os.path.join(input_dataset_dir, "alignments", alignment_name)

  output_species_tree = os.path.join(output, "speciesTree.newick")
  output_true_tree = os.path.join(output, "trueGeneTree.newick")
  output_raxml_tree = os.path.join(output, "raxmlGeneTree.newick")
  output_treerecs_tree = os.path.join(output, "treerecsGeneTree.newick")
  output_msa = os.path.join(output, "alignment.msa")


  os.makedirs(output)
  shutil.copyfile(species_tree, output_species_tree)
  shutil.copyfile(msa_file, output_msa)
  family_index = find_family_index(alignments_info_file, alignment_name)
  print(family_index)
  if (family_index == -1):
    print(alignment_name)
  extract_line(true_trees, family_index, output_true_tree)
  extract_line(raxml_trees, family_index, output_raxml_tree)
  extract_line(treerecs_trees, family_index, output_treerecs_tree)


if __name__ == '__main__':
  if (len(sys.argv) != 4):
    print("Syntax : input_dataset_dir alignment_name output_dir")
    sys.exit(1)

  input_dataset_dir = sys.argv[1]
  alignment_name = sys.argv[2]
  output = sys.argv[3]

  extract_gene_family(input_dataset_dir, alignment_name, output)
