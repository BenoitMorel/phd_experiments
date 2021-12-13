import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import fam
import experiments as exp


def uncompress(archive):
  archive_name = os.path.basename(os.path.normpath(archive)).replace(".tar.gz", "")
  datadir = fam.get_datadir(archive_name)
  families_path = os.path.abspath(fam.get_datasets_family_path())
  if (os.path.isdir(datadir)):
      print("  [ERROR] " + datadir + " already exists. Abording!")
      sys.exit(1)
  print("  [INFO] Starting extracting the archive " + archive)
  os.mkdir(datadir)
  command = []
  command.append("tar")
  command.append("-xzf")
  command.append(archive)
  command.append("-C")
  command.append(families_path)
  subprocess.check_call(command)
  print(" ".join(command))
  if (not os.path.isdir(datadir)):
    print("  [ERROR] Failed to extract the archive to " + datadir)
    sys.exit(1)
  print("  [INFO] Archive successfully extract into " + datadir)

def unarchive(archive):
  if (not os.path.isfile(archive) or not archive.endswith(".tar.gz")):
    print("  [ERROR] " + archive + " is not a valid archive")
    sys.exit(1)
  uncompress(archive)


if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " datadir")
    sys.exit(1)
  unarchive(sys.argv[1])




