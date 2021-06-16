import os
import sys
import shutil
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp
from os.path import expanduser

# sudo apt-get install libgsl0-dev
# setenforce 0
# sudo cpan 
# force install Math::GSL

def call(command):
  print(" ".join(command))
  subprocess.check_call(command)

def mkdir(dir_name):
  try:
    os.makedirs(dir_name)
  except:
    pass

def git_update(repo, output = "", branch = "", output_prefix = exp.github_root):
  cwd = os.getcwd()
  output = os.path.join(output_prefix, output)
  if (os.path.isdir(output)):
    print("Git repository " + output + " already exists")
  else:
    print("Cloning " + repo + " into " + output)
    command = []
    command.append("git")
    command.append("clone")
    command.append("--recursive")
    if (len(branch) > 0):
      command.append("--branch")
      command.append(branch)
    command.append(repo)
    command.append(output)
    call(command)
  os.chdir(output)
  if (len(branch) > 1):
    call(["git", "checkout", branch])
  call(["git", "pull"])
  try:
    call(["git", "submodule", "update", "--init"])
  except:
    print("failed to update submodule for " + repo)
  os.chdir(cwd)

def apply_diff(diff_file, reverse = False):
  cwd = os.getcwd()
  os.chdir(exp.github_root)
  if (reverse):
    os.system("patch -p0 -R < " + diff_file)

  else:
    os.system("patch -p0 < " + diff_file)
  os.chdir(cwd)


def apply_git_diff(repo_name, diff_file):
  output = os.path.join(exp.github_root, repo_name)
  cwd = os.getcwd()
  os.chdir(output)
  call(["git", "checkout", "."])
  call(["git", "apply", os.path.join(exp.installer_root, diff_file)])
  os.chdir(cwd)

def wget(link, file_name, prefix = exp.github_root, unzip = True, unzip_command = ["unzip"]):
  output = os.path.join(prefix, file_name)
  cwd = os.getcwd()
  os.chdir(prefix)
  if (os.path.isfile(output)):
    print("Directory " + output + " already exists")
  else:
    call(["wget", "-O", file_name, link, "-P", exp.github_root])
  unzip_command.append(output)
  if (not os.path.isdir(os.path.splitext(output)[0]) and unzip):
    call(unzip_command)
  os.chdir(cwd)

def install_simphy():
  wget("https://github.com/adamallo/SimPhy/releases/download/v1.0.2/SimPhy_1.0.2.tar.gz", "SimPhy.tar.gz", unzip = True, unzip_command =["tar", "-xzf"])
  call(["chmod", "777", "../SimPhy_1.0.2/bin/simphy_lnx64"])

def install_with_cmake(repo_name, cmake_additional_commands = [], install = False):
  output = os.path.join(exp.github_root, repo_name)
  cwd = os.getcwd()
  os.chdir(output)
  mkdir("build")
  os.chdir("build")
  cmake_command = ["cmake"]
  cmake_command.extend(cmake_additional_commands)
  cmake_command.append("..")
  call(cmake_command)
  call(["cmake", ".."])
  try:
    call(["make", "-j", "40"])
  except:
    call(["make"])

  if (install):
    call(["make", "install"])
  os.chdir(cwd)

def install_with_install_script(repo_name):
  output = os.path.join(exp.github_root, repo_name)
  cwd = os.getcwd()
  os.chdir(output)
  call(["/bin/bash", "install.sh"])
  os.chdir(cwd)

def install_seq_gen(repo):
  output = os.path.join(exp.github_root, repo_name)
  cwd = os.getcwd()
  os.chdir(os.path.join(output, "source"))
  call(["make"])
  os.chdir(cwd)



  





###################
###   MAIN      ###
###################
if (True):
  git_update("https://github.com/BenoitMorel/GeneRax.git", "GeneRax")
  git_update("https://github.com/BenoitMorel/Pargenes.git", "pargenes")
  git_update("https://github.com/amkozlov/raxml-ng.git", "raxml-ng", "dev")
  install_with_cmake("GeneRax")
  install_with_install_script("pargenes") 
  git_update("https://github.com/BenoitMorel/MPIScheduler.git", "MPIScheduler")
  install_with_cmake("MPIScheduler")
  wget("https://github.com/rambaut/Seq-Gen/archive/1.3.4.zip", "seq-gen.zip")
  install_seq_gen("Seq-Gen-1.3.4") 



