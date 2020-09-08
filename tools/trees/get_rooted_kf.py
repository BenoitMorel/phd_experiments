import sys
import os
sys.path.insert(0, 'scripts')
import experiments as exp
import subprocess

  

def get_rooted_kf(tree1, tree2):
    command = []
    command.append("R")
    command.append("-f")
    command.append(exp.getrootedkf_R_script)
    command.append("--args")
    command.append(tree1)
    command.append(tree1)
    subprocess.check_call(command)


if __name__ == '__main__':
    if (len(sys.argv) != 3):
        print("Syntax:  tree1 tree1")
        sys.exit(1)
    tree1 = sys.argv[1]
    tree2 = sys.argv[2]
    get_rooted_kf(tree1, tree1)
    
