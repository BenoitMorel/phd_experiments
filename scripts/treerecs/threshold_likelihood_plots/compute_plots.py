import os
import subprocess

reporoot = "./"
treerecsroot = os.path.join("..", "Treerecs")
datadir = os.path.join(reporoot, "datasets", "phyldog_test_dataset", "treerecs_format")
gene = "HBG011000"
genedir = os.path.join(datadir, gene)
resultsdir= os.path.join(reporoot, "results", "treerecs", "threshold_likelihood_plots", gene)
os.makedirs(resultsdir, exist_ok=True)
script_path = os.path.join(treerecsroot, "scripts", "experiments", "plot_likelihood_vs_threshold.py")

command = []
command.append("python")
command.append(script_path)
command.append(os.path.join(genedir, gene + "_raxml_support.newick"))
command.append(os.path.join(datadir, "speciesTree.newick"))
command.append(os.path.join(genedir, "alignment.txt"))
command.append(os.path.join(genedir, gene + ".map"))
command.append("7")
command.append(resultsdir)
subprocess.check_call(command)
#python $treerecsroot/scripts/experiments/plot_likelihood_vs_threshold.py ${datadir}/${gene}/${gene}_raxml_support.newick ${datadir}/speciesTree.newick ${datadir}/${gene}/alignment.txt ${datadir}/${gene}/${gene}.map 20 $resultsdir

