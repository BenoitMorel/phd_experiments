import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import fam
import experiments as exp


def compress(datadir):
  datadir = os.path.normpath(datadir)
  dataset = os.path.basename(datadir)
  output = datadir + ".tar.gz"
  if (os.path.isfile(output)):
    print("  [INFO] Removing existing archive " + output)
    os.remove(output)
  print("  [INFO] Starting creating archive " + output)
  families_path = os.path.abspath(fam.get_datasets_family_path())
  command = []
  command.append("tar")
  command.append("-czf")
  command.append(output)
  command.append("--directory=" + families_path)
  command.append(dataset)
  print(" ".join(command))
  subprocess.check_call(command)
  if (not os.path.isfile(output)):
    print("  [ERROR] the archive " + output + " was not produced...")
    sys.exit(1)
  print("  [INFO] Archive " + output + " successfully created")
  return output

def backup(output):
  output_name = os.path.basename(os.path.normpath(output))
  backup_path = os.path.join(exp.fast_dataset_archive, output_name)
  print("  [INFO[ Creating backup in " + backup_path)
  shutil.copyfile(output, backup_path)

def clean(datadir):
  print("  [INFO] Removing the initial directory " + datadir)
  shutil.rmtree(datadir)

def archive(datadir):
  if (not os.path.isdir(datadir)):
    print("  [ERROR] " + datadir + " is not an archivable directory")
    sys.exit(1)
  output = compress(datadir)
  backup(output)
  clean(datadir)


if (__name__ == "__main__"):
  if (len(sys.argv) < 2):
    print("Syntax python " + os.path.basename(__file__) + " datadir_list")
    sys.exit(1)
  datadir_list = sys.argv[1:]
  for datadir in datadir_list:
    archive(datadir)



