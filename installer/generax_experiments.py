import os
import sys
import shutil
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp

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

def download_jprime(repo_name):
  output = os.path.join(exp.github_root, repo_name)
  cwd = os.getcwd()
  call(["wget", "-O", "jprime-0.3.6.jar", "https://uc576760c336faa30a378bf489e0.dl.dropboxusercontent.com/cd/0/get/AeqjqdemOXeK2NCiK5V61mMDRn-wt-Pt_8moifOti9pcBQHRZbhUfA6pu8ihBiluEjrlAYt5yS5T3XuOoVfdP1tJmCl8W7cPlyNwuixWu8SpIQ/file?_download_id=42561867786005014967751414617589329163419213240565884563046145126"])
  os.chdir(cwd)

def install_seq_gen(repo):
  output = os.path.join(exp.github_root, repo_name)
  cwd = os.getcwd()
  os.chdir(os.path.join(output, "source"))
  call(["make"])
  os.chdir(cwd)


def install_bpp_for_ale():
  call([os.path.join(exp.installer_root, "data", "bpp-setup.sh")])

def add_string_to_file(f, s):
  content = open(f).read()
  with open(f, "w") as writer:
    writer.write(content + "\n" + s + "\n")

def install_standard_raxml(name):
  cwd = os.getcwd()
  os.chdir(name) 
  call(["make", "Makefile.AVX.gcc"])
  os.chdir(cwd)

def install_phyldog(repo_name):
  cwd = os.getcwd()
  phyldog_root = os.path.join(exp.github_root, repo_name)
  phyldog_deps = os.path.join(phyldog_root, "temp_install")
  #mkdir(phyldog_deps)
  boost_lib_path = os.path.join(phyldog_deps, "boost_1_58_0", "boost_install", "lib")
  boost_include_path = os.path.join(phyldog_deps, "boost_1_58_0", "boost_install", "include")
  bpp_install = os.path.join(phyldog_deps, "bpp_install")
  bpp_libs = os.path.join(bpp_install, "lib")
  bpp_include = os.path.join(bpp_install, "include")
  #apply_git_diff("PHYLDOG", "phyldog_diff.txt")
  os.chdir(phyldog_deps)
  #wget("https://sourceforge.net/projects/boost/files/boost/1.58.0/boost_1_58_0.zip/download", "boost.zip", phyldog_deps)
  os.chdir("boost_1_58_0")
  #call(["./bootstrap.sh", "--prefix=boost_install", "--with-libraries=mpi,serialization"])
  #add_string_to_file("project-config.jam", "using mpi ;")
  #call(["./b2", "--user-config=project-config.jam"]) 
  #call(["./b2", "install"]) 
  
 
  os.chdir(os.path.join(phyldog_deps))
  #git_update("https://github.com/BioPP/bpp-core.git", "bpp-core", "rel220", ".")
  #git_update("https://github.com/BioPP/bpp-seq.git", "bpp-seq", "rel220", ".")
  #git_update("https://github.com/BioPP/bpp-phyl.git", "bpp-phyl", "rel220", ".")

  #install_with_cmake(os.path.join(phyldog_deps, "bpp-core"), ["-DCMAKE_INSTALL_PREFIX=" + bpp_install], True)
  #install_with_cmake(os.path.join(phyldog_deps, "bpp-seq"), ["-DCMAKE_INSTALL_PREFIX=" + bpp_install], True)
  #install_with_cmake(os.path.join(phyldog_deps, "bpp-phyl"), ["-DCMAKE_INSTALL_PREFIX=" + bpp_install, "-DCMAKE_LIBRARY_PATH=" + bpp_libs, "-DCMAKE_INCLUDE_PATH=" + bpp_include], True)
  print ("WARNING: YOU NEED TO DOWNLOAD AND BUILD OLD PLL HERE!!!")
  pll_libs = os.path.join(phyldog_deps, "pll") 
  pll_includes = os.path.join(phyldog_deps, "pll") 
  os.chdir(phyldog_root)
  
  all_libs = boost_lib_path + ";" + bpp_libs + ";" + pll_libs
  all_includes = boost_include_path + ";" + bpp_include + ";" + pll_includes
  install_with_cmake(phyldog_root, ["-DCMAKE_LIBRARY_PATH=" + all_libs, "-DCMAKE_INCLUDE_PATH=" + all_includes]) #, "-DBUILD_STATIC=ON"])
  
  os.chdir(cwd)
  
  

def run_make(repo):
  os.chdir(repo)
  call(["make"])
  os.chdir(cwd)

def install_deco(targz):
  cwd = os.getcwd()
  os.chdir(exp.github_root)
  subprocess.check_call(["gunzip", targz])
  subprocess.check_call(["tar", "-xf", targz[:-3]])
  repo = os.path.join(exp.github_root, targz)[:-7]
  print(repo)
  
  apply_diff(os.path.join(exp.github_root, "phd_experiments", "installer", "deco_diff.txt"), reverse = True)
  apply_diff(os.path.join(exp.github_root, "phd_experiments", "installer", "deco_make_diff.txt"), reverse = True)
  run_make(repo)

if (False):
  git_update("https://github.com/ssolo/ALE.git", "ALE")
  git_update("https://github.com/BenoitMorel/BenoitDatasets.git", "BenoitDatasets")
  git_update("https://github.com/aberer/exabayes.git", "exabayes-1.5")
  git_update("https://github.com/BenoitMorel/GeneRax.git", "GeneRax")
  download_jprime("jprime")
  wget("http://goby.compbio.cs.cmu.edu/Notung/Notung-2.9.zip", "Notung-2.9.zip")
  git_update("https://github.com/BenoitMorel/Pargenes.git", "pargenes")
  git_update("https://github.com/amkozlov/raxml-ng.git", "raxml-ng", "dev")
  git_update("https://gitlab.inria.fr/Phylophile/Treerecs.git", "Treerecs", "treesearch")
  git_update("https://github.com/stamatak/standard-RAxML.git", "standard-RAxML")  
  install_bpp_for_ale()
  apply_git_diff("ALE", "ale_diff.txt")
  home = os.path.expanduser("~")
  install_with_cmake("ALE", ["-DCMAKE_LIBRARY_PATH=" + home + "/install/bio++/lib", "-DCMAKE_INCLUDE_PATH=" + home + "/install/bio++/include/"])
  
  wget("https://cme.h-its.org/exelixis/resource/download/software/exabayes-1.5.zip", "exabayes-1.5.zip")
  apply_diff(os.path.join(exp.github_root, "phd_experiments", "installer", "exabayes_diff.txt"))
  install_with_cmake("GeneRax")
  install_with_install_script("pargenes") 
  git_update("https://github.com/BenoitMorel/MPIScheduler.git", "MPIScheduler")
  install_with_cmake("MPIScheduler")
  install_with_cmake("raxml-ng")
  install_with_cmake("Treerecs")
  wget("https://github.com/rambaut/Seq-Gen/archive/1.3.4.zip", "seq-gen.zip")
  install_seq_gen("Seq-Gen-1.3.4") 

  treerecs_exec = os.path.join(exp.github_root, "Treerecs", "build", "bin", "treerecs")
  treerecs_exec_2 = os.path.join(exp.github_root, "Treerecs", "build", "bin", "treerecs")
  try:
    shutil.copy(treerecs_exec, treerecs_exec_2)
  except:
    pass
  install_standard_raxml("standard-RAxML")

  git_update("https://github.com/Boussau/PHYLDOG", "PHYLDOG")

  install_phyldog("PHYLDOG")
  install_with_autotools("exabayes-1.5")
  wget("http://pbil.univ-lyon1.fr/software/DeCo/DeCo.tar.gz", "DeCo.tar.gz", unzip = False)
  install_deco("DeCo.tar.gz")
  
  git_update("https://github.com/WandrilleD/DeCoSTAR.git", "DeCoSTAR")
  subprocess.check_call(["./installer/install_recent_bpp.sh"], shell = True)
  apply_diff(os.path.join(exp.github_root, "phd_experiments", "installer", "decostart_make_diff.txt"))#, reverse = True)
  git_update("https://github.com/davidemms/STAG.git", "STAG")
  
  git_update("https://github.com/celinescornavacca/ecceTERA.git", "ecceTERA")
  install_with_cmake("ecceTERA")

if (True):
  #git_update("https://github.com/Boussau/PHYLDOG", "PHYLDOG")
  install_phyldog("PHYLDOG")
  #install_simphy()


