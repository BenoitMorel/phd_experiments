import os
import sys
import shutil
import subprocess





class Common():
  def __init__(self, git_dir, cores):
    self.git_dir = git_dir
    self.installer_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
    self.cores = cores


  def call(self, command):
    print(" ".join(command))
    subprocess.check_call(command)

  def mkdir(self, dir_name):
    try:
      os.makedirs(dir_name)
    except:
      pass

  def git_update(self, repo, output = "", branch = "", output_prefix = None):
    if (output_prefix == None):
      output_prefix = self.git_dir
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
      self.call(command)
    os.chdir(output)
    if (len(branch) > 1):
      self.call(["git", "checkout", branch])
    self.call(["git", "pull"])
    try:
      self.call(["git", "submodule", "update", "--init"])
    except:
      print("failed to update submodule for " + repo)
    os.chdir(cwd)

  def apply_diff(self, diff_file, reverse = False):
    cwd = os.getcwd()
    os.chdir(self.git_dir)
    if (reverse):
      os.system("patch -p0 -R < " + diff_file)
    else:
      os.system("patch -p0 < " + diff_file)
    os.chdir(cwd)

  def apply_git_diff(self, repo_name, diff_file):
    output = os.path.join(self.git_dir, repo_name)
    cwd = os.getcwd()
    os.chdir(output)
    self.call(["git", "checkout", "."])
    self.call(["git", "apply", os.path.join(self.installer_path, diff_file)])
    os.chdir(cwd)

  def wget(self, link, file_name, prefix = None, unzip = True, unzip_command = ["unzip"]):
    if (prefix == None):
      prefix = self.git_dir
    output = os.path.join(prefix, file_name)
    cwd = os.getcwd()
    os.chdir(prefix)
    if (os.path.isfile(output)):
      print("Directory " + output + " already exists")
    else:
      self.call(["wget", "-O", file_name, link, "-P", self.git_dir])
    unzip_command.append(output)
    if (not os.path.isdir(os.path.splitext(output)[0]) and unzip):
      self.call(unzip_command)
    os.chdir(cwd)

  def install_with_cmake(self, repo_name, cmake_additional_commands = [], install = False):
    output = os.path.join(self.git_dir, repo_name)
    cwd = os.getcwd()
    os.chdir(output)
    self.mkdir("build")
    os.chdir("build")
    cmake_command = ["cmake"]
    cmake_command.extend(cmake_additional_commands)
    cmake_command.append("..")
    self.call(cmake_command)
    self.call(["cmake", ".."])
    try:
      self.call(["make", "-j", str(self.cores)])
    except:
      self.call(["make"])

    if (install):
      self.call(["make", "install"])
    os.chdir(cwd)

  def install_with_autotools(self, repo_name):
    output = os.path.join(self.git_dir, repo_name)
    cwd = os.getcwd()
    os.chdir(output)
#call(["./bootstrap.sh"])
    self.call(["./configure"])
    self.call(["make", "-j", str(self.cores)])
    os.chdir(cwd)


  def get_installer_path(self):
    return self.installer_path




