import os
import sys
import subprocess
import shutil
import time
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "msa_edition"))
import fam
import saved_metrics
import experiments as exp
import msa_converter
import rename_leaves
import sequence_model
import run_raxml_supportvalues as run_pargenes
import read_tree




def generate_config_file(chains, generations, frequency, output_config_file, nexus_alignment, subst_model, seed, output_prefix):
    with open(output_config_file, "w") as writer:
      writer.write("\tbegin mrbayes;\n")
      writer.write("\tset seed=" + str(seed) + ";\n")
      writer.write("\tset autoclose=yes nowarn=yes;\n")
      writer.write("\texecute " + nexus_alignment + ";\n")
      writer.write(sequence_model.get_mrbayes_preset_line(subst_model))
      writer.write(sequence_model.get_mrbayes_lset_line(subst_model))
      writer.write("\tmcmc nruns=1" + " nchains=" + str(chains) + " ngen=" + str(generations) + " samplefreq=" + str(frequency) + " file=" + output_prefix + " append=no" + ";\n")
      writer.write("end;")


def generate(alignment, subst_model, runs, chains, generations, frequency, burnin, outputdir):
    os.mkdir(outputdir)
    prefix = os.path.join(outputdir, "mrbayes")
    config_file = os.path.join(outputdir, "config_file.txt")
    # create nexus alignment
    nexus_alignment = os.path.join(outputdir, "alignment.nexus")
    msa_converter.msa_convert(alignment, nexus_alignment, "fasta", "nexus", None)
    # create mrbayes config file
    generate_config_file(chains, generations, frequency, config_file, nexus_alignment, subst_model, 42, prefix) 
    # generate the command
    command = []
    command.append(exp.mrbayes_exec)
    command.append(config_file)
    subprocess.run(command)


if (__name__== "__main__"):
  if len(sys.argv) != 9:
    print("Syntax error: python " + os.path.basename(__file__) + " alignment subst_model runs chains generations frequency burnin output.newick")
    print(len(sys.argv))
    sys.exit(0)

  alignment = sys.argv[1]
  subst_model = sys.argv[2]
  runs = int(sys.argv[3])
  chains = int(sys.argv[4])
  generations = int(sys.argv[5])
  frequency = int(sys.argv[6])
  burnin = sys.argv[7]
  output = sys.argv[8]
  generate(alignment, subst_model, runs, chains, generations, frequency, burnin, output)
  

