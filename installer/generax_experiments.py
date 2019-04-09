import os
import sys
import shutil
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp



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

def apply_git_diff(repo_name, diff_file):
  output = os.path.join(exp.github_root, repo_name)
  cwd = os.getcwd()
  os.chdir(output)
  call(["git", "checkout", "."])
  call(["git", "apply", os.path.join(exp.installer_root, diff_file)])
  os.chdir(cwd)

def wget(link, file_name):
  output = os.path.join(exp.github_root, file_name)
  cwd = os.getcwd()
  os.chdir(exp.github_root,)
  if (os.path.isfile(output)):
    print("Directory " + output + " already exists")
  else:
    call(["wget", link, "-P", exp.github_root])
  if (not os.path.isdir(os.path.splitext(output)[0])):
    call(["unzip", output])
  os.chdir(cwd)

def install_with_cmake(repo_name, cmake_additional_commands = []):
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
  call(["make", "-j", "40"])
  os.chdir(cwd)

def install_with_autotools(repo_name):
  output = os.path.join(exp.github_root, repo_name)
  cwd = os.getcwd()
  os.chdir(output)
  #call(["./bootstrap.sh"])
  call(["./configure"])
  call(["make", "-j", "40"])
  os.chdir(cwd)

def install_with_install_script(repo_name):
  output = os.path.join(exp.github_root, repo_name)
  cwd = os.getcwd()
  os.chdir(output)
  call(["/bin/bash", "install.sh"])
  os.chdir(cwd)

def install_bpp_for_ale():
  
  call([os.path.join(exp.installer_root, "data", "bpp-setup.sh")])
  
if (False):
  git_update("https://github.com/ssolo/ALE.git", "ALE")
  git_update("https://github.com/BenoitMorel/BenoitDatasets.git", "BenoitDatasets")
  git_update("https://github.com/aberer/exabayes.git", "exabayes-1.5")
  git_update("https://github.com/BenoitMorel/GeneRax.git", "GeneRax")
  git_update("https://github.com/arvestad/jprime.git", "jprime")
  wget("http://goby.compbio.cs.cmu.edu/Notung/Notung-2.9.zip", "Notung-2.9.zip")
  git_update("https://github.com/BenoitMorel/Pargenes.git", "pargenes")
  git_update("https://github.com/amkozlov/raxml-ng.git", "raxml-ng", "dev")
  git_update("https://gitlab.inria.fr/Phylophile/Treerecs.git", "Treerecs", "dev")
  
  install_bpp_for_ale()
  apply_git_diff("ALE", "ale_diff.txt")
  home = os.path.expanduser("~")
  install_with_cmake("ALE", ["-DCMAKE_LIBRARY_PATH=" + home + "/install/bio++/lib", "-DCMAKE_INCLUDE_PATH=" + home + "/install/bio++/include/"])
  
  wget("https://cme.h-its.org/exelixis/resource/download/software/exabayes-1.5.zip", "exabayes-1.5.zip")
  install_with_autotools("exabayes-1.5")
  install_with_cmake("GeneRax")
  install_with_install_script("pargenes") 
  git_update("https://github.com/BenoitMorel/MPIScheduler.git", "MPIScheduler")
  install_with_cmake("MPIScheduler")
  install_with_cmake("raxml-ng")

if (True):
  
  install_with_cmake("Treerecs")
  
  
  #git_update("", "")
  #git_update("", "")
  #git_update("", "")
  #git_update("", "")
  #git_update("", "")






