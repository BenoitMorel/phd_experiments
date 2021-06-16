import os
import common
from common import Common

class Installer():
  def __init__(self, git_dir, cores):
    self.common = Common(git_dir, cores)
    self.common.mkdir(git_dir)
    self.git_dir = git_dir
    self.cores = cores
    
  def install_astralpro(self):
    self.common.git_update("https://github.com/chaoszhang/A-pro.git", "A-pro")
  
  def install_generax(self):
    self.common.git_update("https://github.com/BenoitMorel/GeneRax.git", "GeneRax", "genetegrator")
    self.common.install_with_cmake("GeneRax")

  def install_mpischeduler(self):
    self.common.git_update("https://github.com/BenoitMorel/MPIScheduler.git", "MPIScheduler")
    self.common.install_with_cmake("MPIScheduler")
  
  def install_mrbayes(self):
    self.common.git_update("https://github.com/NBISweden/MrBayes.git", "MrBayes")
    self.common.install_with_autotools("MrBayes")
  
  def install_raxmlng(self):
    self.common.git_update("https://github.com/amkozlov/raxml-ng.git", "raxml-ng", "dev")
    self.common.install_with_cmake("raxml-ng")

  def install_seq_gen(self):
    self.common.wget("https://github.com/rambaut/Seq-Gen/archive/1.3.4.zip", "seq-gen.zip")
    output = os.path.join(self.git_dir, "Seq-Gen-1.3.4")
    cwd = os.getcwd()
    os.chdir(os.path.join(output, "source"))
    self.common.call(["make"])
    os.chdir(cwd)

  def install_simphy(self):
    installer_path = self.common.get_installer_path()
    simphy_installer = os.path.join(installer_path, "install_simphy.sh")
    self.common.call([simphy_installer, self.git_dir])

