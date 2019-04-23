import sys
import os
import subprocess
import re

def query_haswell_vpn(user_name, query):
  command = []
  command.append("ssh")
  command.append("-t")
  command.append(user_name + "@bridge-login.h-its.org")
  command2 = []
  command2.append("ssh")
  command2.append("-t")
  command2.append(user_name + "@cluster-login-hits.urz.uni-heidelberg.de")
  command2.append(query)
  command.append(" ".join(command2))
  return subprocess.check_output(command)

def query_haswell(user_name, query):
  command = []
  command.append("ssh")
  command.append("-t")
  command.append(user_name + "@cluster-login-hits.urz.uni-heidelberg.de")
  command.append(query)
  return subprocess.check_output(command)
  

def get_running_jobs(user_name):
  logs = query_haswell_vpn(user_name, "squeue")
  jobs = {}
  for line in logs.split("\n"):
    if (("haswell" in line) and (user_name in line)):
      line = re.sub(' +', ' ', line)
      jobs[line.split(" ")[1]] = True
  return jobs

def store_jobs(jobs):
  temp_file = os.path.join(os.path.expanduser("~"), ".temp_haswell_notifier.txt")
  with open(os.path.join(temp_file), "w") as writer:
    for job in jobs:
      writer.write(job + " ")

def get_last_read_jobs():
  temp_file = os.path.join(os.path.expanduser("~"), ".temp_haswell_notifier.txt")
  jobs = {}
  try:
    line = open(os.path.join(temp_file)).read()
    for job in line.split(" "):
      if (len(job)):
        jobs[job] = True
  except:
    pass
  return jobs


previous_jobs = get_last_read_jobs()
print(" last jobs: " + str(previous_jobs))
current_jobs = get_running_jobs("morelbt")
print(" current jobs: " + str(current_jobs))
store_jobs(current_jobs)
