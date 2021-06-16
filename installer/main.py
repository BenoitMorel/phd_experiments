import installer
import sys
import os

if (__name__ == "__main__"):
  if (len(sys.argv) < 2):  
    print("Syntax python " + os.path.basename(__file__) + "install_dir [cores=2]")
    sys.exit(1)
  git_dir = sys.argv[1]
  cores = 2
  if (len(sys.argv) > 2):
    cores = int(sys.argv[2])

  inst = installer.Installer(git_dir, cores)
  inst.install_raxmlng()
  inst.install_seq_gen()
  inst.install_simphy()
  inst.install_astralpro()
  inst.install_generax()
  inst.install_mrbayes()
  inst.install_mpischeduler()

