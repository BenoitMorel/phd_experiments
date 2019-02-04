import sys
import os
import shutil
sys.path.insert(0, 'tools/msa_edition')
sys.path.insert(0, 'tools/raxml')
import msa_subsampler
import raxml_get_tca_score


def get_output_dir(input_dataset, trees_selection_ratio, sampling_ratio):
  parent_dir = os.path.abspath(os.path.join(input_dataset, os.pardir))
  output_dataset = "sub_t" + str(trees_selection_ratio) + "_s" + str(sampling_ratio) + "_" 
  output_dataset += os.path.basename(os.path.normpath(input_dataset))
  output_dataset = os.path.join(parent_dir, output_dataset)
  return output_dataset

def get_selected_families(input_families_dir, trees_selection_ratio):
  family_tcas = {}
  for family in os.listdir(input_families_dir):
    try:
      tca = float(open(os.path.join(input_families_dir, family, "tca.txt")).readline())
    except:
      continue
    family_tcas[family] = tca
  sorted_families = sorted(family_tcas, key=family_tcas.get, reverse=True)
  final_size = int(float(len(sorted_families)) * trees_selection_ratio)
  if (final_size == 0):
    print("Error: 0 families were selected")
    exit(1)
  return sorted_families[0:final_size]


def extract_main_files(input_dataset, output_dataset):
  input_species = os.path.join(input_dataset, "speciesTree.newick")
  output_species = os.path.join(output_dataset, "speciesTree.newick")
  shutil.copyfile(input_species, output_species)

def mycopy(dir1, dir2, file1, file2 = None):
  if (file2 == None):
    file2 = file1
  shutil.copy(os.path.join(dir1, file1), os.path.join(dir2, file2))


def extract_family(input_family, output_family, output_alignments_dir, sampling_ratio):
  os.makedirs(output_family)
  family = os.path.basename(os.path.normpath(input_family))
  mycopy(input_family, output_family, "mapping.link") 
  mycopy(input_family, output_family, "treerecs_mapping.link") 
  mycopy(input_family, output_family, "raxmlGeneTree.newick", "trueGeneTree.newick") 
  input_alignment = os.path.join(input_family, "alignment.msa")
  output_alignment = os.path.join(output_family, "alignment.msa")
  msa_format = "fasta"
  msa_subsampler.subsample(input_alignment,output_alignment, msa_format, sampling_ratio)
  shutil.copy(output_alignment, os.path.join(output_alignments_dir, family + "." + msa_format))




  treerec_dir = os.path.join(output_family, "treerecs")
  os.makedirs(treerec_dir)
  treerecs_model = open(os.path.join(input_family, "treerecs", "alignment_descriptor.txt")).readline()
  with open(os.path.join(treerec_dir, "alignment_descriptor.txt"), "w") as writer:
    writer.write(treerecs_model)
    writer.write(os.path.abspath(output_alignment))

def subsample_best_tree(input_dataset, trees_selection_ratio, sampling_ratio):
  # selection
  input_families_dir = os.path.join(input_dataset, "families")
  selected_families = get_selected_families(input_families_dir, trees_selection_ratio)
  
  # start extracting selected families
  output_dataset = get_output_dir(input_dataset, trees_selection_ratio, sampling_ratio)
  os.makedirs(output_dataset)
  output_families_dir = os.path.join(output_dataset, "families")
  os.makedirs(output_families_dir)
  extract_main_files(input_dataset, output_dataset)
  output_alignments_dir = os.path.join(output_dataset, "alignments")
  os.makedirs(output_alignments_dir)
  for family in selected_families:
    input_family_dir = os.path.join(input_families_dir, family)
    output_family_dir = os.path.join(output_families_dir, family)
    extract_family(input_family_dir, output_family_dir, output_alignments_dir, sampling_ratio)
    mycopy(input_dataset, output_family_dir, "speciesTree.newick")
  print("Output directory: " + output_dataset)


if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax: python subsample_best_trees.py input_dataset  trees_selection_ratio sampling_ratio")
    exit(1)
  input_dataset = sys.argv[1]
  trees_selection_ratio = float(sys.argv[2])
  sampling_ratio = float(sys.argv[3])
  subsample_best_tree(input_dataset, trees_selection_ratio, sampling_ratio)





