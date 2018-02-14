import sys
import os
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.insert(0, 'scripts')
import experiments as exp

def run_treerecs(gene_tree, species_tree, alignment, smap, thresholds_number, output_dir):
    treerecs_exec = exp.treerecs_exec
    treerecs_output = os.path.join(output_dir, "treerecs_output.newick")
    command = []
    command.append(treerecs_exec)
    command.append("-g")
    command.append(gene_tree)
    command.append("-s")
    command.append(species_tree)
    command.append("-o")
    command.append(treerecs_output)
    command.append("-a")
    command.append(alignment)
    command.append("-t")
    command.append("all")
    command.append("--ale-evaluation")
    command.append("-T")
    command.append(thresholds_number)
    command.append("-S")
    command.append(smap)
    print(' '.join(command))
    subprocess.check_call(command)
    return treerecs_output

#
# ll_dico[family][threshold][ll_name] = ll
#   
# t[family]
# ll[family[llname]
def fill_dico(lines, thresholds, likelihoods):
    for line in lines:
        if (not line.startswith(">")):
            continue
        split = line.split(" ")
        family_name = ""
        threshold = ""
        
        entry = {}
        for index, val in enumerate(split):
            if (val == "family"):
                family_name = split[index + 1]
                if (thresholds.get(family_name) == None):
                    thresholds[family_name] = []
                    likelihoods[family_name]= {}
            elif (val == "contraction" and split[index + 1] == "threshold"):
                threshold = split[index + 3][:-1]
                if (threshold == "no"):
                    threshold = "0.0"
                threshold = float(threshold)
            elif (val == "logLk"):
                entry[split[index - 1]] = float(split[index + 2][:-1])
        
        thresholds[family_name].append(threshold)
        join_ll = 0
        for llname, llvalue in entry.items():
            if (likelihoods[family_name].get(llname) == None):
                likelihoods[family_name][llname] = []
            likelihoods[family_name][llname].append(llvalue)
            join_ll += float(llvalue)
        if (likelihoods[family_name].get("joint") == None):
            likelihoods[family_name]["joint"] = []
        likelihoods[family_name]["joint"].append(join_ll)

        
def export_family(family_name, thresholds, likelihoods, output_dir):
    sns.set()
    output = os.path.join(output_dir, "family_" + family_name)
    fig = plt.figure()
    index = 311
    for llname, llvalues in likelihoods[family_name].items():
        plt.subplot(index)
        index = index + 1
        plt.plot(thresholds[family_name], llvalues)
        legend = []
        legend.append(llname)
        plt.xlabel('Thresholds')
        plt.ylabel('LogLikelihood')
        plt.legend(legend, loc='upper left')
        plt.tight_layout()
       
    fig.savefig(output)

def parse_treerecs_output(treerecs_output, output_dir):
    thresholds = {}
    likelihoods = {}
    with open(treerecs_output) as f:
        fill_dico(f.readlines(), thresholds, likelihoods)
    for family in thresholds:
        export_family(family, thresholds, likelihoods, output_dir)



def compute_plots(datadir, resultsdir):
  gene_trees = os.path.join(datadir, "geneTrees.newick")
  species_trees = os.path.join(datadir, "speciesTree.newick")
  alignment = os.path.join(datadir, "alignment.txt")
  smap = os.path.join(datadir, "mapping.txt")
  thresholds_number = "7"
  run_treerecs(gene_trees, species_trees, alignment, smap, thresholds_number, resultsdir)
  treerecs_output = os.path.join(resultsdir, "treerecs_output.newick")
  parse_treerecs_output(treerecs_output, resultsdir)





datadir = os.path.join(exp.datasets_root, "treerecs")

if (len(sys.argv) != 2):
  datasets = os.listdir(datadir)
  print("Syntax error: dataset required. Suggestions of datasets: ")
  print('\n'.join(datasets))
  sys.exit(0)

basedir = sys.argv[1]
datadir = os.path.join(datadir, basedir)
resultsdir = exp.create_result_dir(os.path.join("treerecs", "threshold_likelihoods_plots", basedir))
result_msg = "Treerecs git: \n" + exp.get_git_info(exp.treerecs_root)
exp.write_results_info(resultsdir, result_msg) 
compute_plots(datadir, resultsdir)





