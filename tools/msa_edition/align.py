import os
import sys
sys.path.insert(0, 'scripts')
import experiments as exp
import subprocess

def align(input_msa, output_msa, aligner):
  command = []
  command.append(exp.mafft_exec)
  command.append(input_msa)
  print(input_msa)
  with open(output_msa, "w") as writer:
    subprocess.check_call(command, stdout=writer)


if (__name__== "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax: input_msa output_msa alignmer")
    print("aligner can be {MAFFT}")
    exit(1)
  align(sys.argv[1], sys.argv[2], sys.argv[3])

