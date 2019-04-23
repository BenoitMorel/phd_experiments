import sys
import os
import subprocess
import re


def get_running_jobs(user_name):
  command = []
  command.append("ssh")
  command.append("-t")
  command.append(user_name + "@cluster-login-hits.urz.uni-heidelberg.de")
  command.append("squeue")
  logs = subprocess.check_output(command)
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
  line = open(os.path.join(temp_file)).read()
  jobs = {}
  for job in line.split(" "):
    if (len(job)):
      jobs[job] = True
  return jobs


previous_jobs = get_last_read_jobs()
print(" last jobs: " + str(previous_jobs))
current_jobs = get_running_jobs("morelbt")
print(" current jobs: " + str(current_jobs))
store_jobs(current_jobs)
