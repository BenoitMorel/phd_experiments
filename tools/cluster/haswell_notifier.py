import sys
import os
import subprocess
import re
import pickle
import time 
import datetime as dt

is_remote = False # set to True if you need to access haswell haswell through magny
user_name = "morelbt" # haswell username
sound_on = True
timer_interval = 10 # elapsed time between queries (s)

def query_haswell_remote(user_name, query):
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
  FNULL = open(os.devnull, 'w')
  return subprocess.check_output(command, stderr=FNULL)

def query_haswell(user_name, query, is_remote):
  if (is_remote):
    return query_haswell_remote(user_name, query)
  command = []
  command.append("ssh")
  command.append("-t")
  command.append(user_name + "@cluster-login-hits.urz.uni-heidelberg.de")
  command.append(query)
  FNULL = open(os.devnull, 'w')
  return subprocess.check_output(command, stderr=FNULL)
  

def get_running_jobs(user_name, is_remote):
  logs = query_haswell(user_name, "squeue", is_remote)
  jobs = {}
  for line in logs.split("\n"):
    if (("haswell" in line) and (user_name in line)):
      line = re.sub(' +', ' ', line)
      split = line.split(" ")
      job_id = split[1]
      status = split[8]
      jobs[job_id] = status
  return jobs


def store_jobs(jobs):
  temp_file = os.path.join(os.path.expanduser("~"), ".temp_haswell_notifier.txt")
  pickle.dump(jobs, open(temp_file, "wb"))

def get_last_read_jobs():
  temp_file = os.path.join(os.path.expanduser("~"), ".temp_haswell_notifier.txt")
  return pickle.load(open(temp_file, "rb"))


def notify(message):
  if (len(message)):
    message.replace("\n", "\\\\n")
    subprocess.Popen(['notify-send', "Haswell-notif",  message], stdin=None, stdout=None, stderr=None, close_fds=True)
  
def play_sounds(sounds):
  try:
    if (sound_on):
      for sound in sounds:
        subprocess.Popen(['spd-say', sound])
        time.sleep(1) #otherwise the sounds overlapp
  except:
    pass


lt_announced = False
def check_lt():
  global lt_announced
  t = dt.datetime.now()
  if (not lt_announced):
    lt_announced = True
    if (t.hour == 12 and t.minute == 25):
      msg = "Warning! Warning! It's lunchtime in 5 minutes!"
      notify(msg)
      play_sounds([msg])
    elif (t.hour == 12 and t.minute == 30):
      msg = "It's lunchtime! Go go go!"
      notify(msg)
      play_sounds([msg])
  else:
    lt_announced = False

def compare_jobs(previous_jobs, new_jobs):
  notifications = []
  sounds = []
  for job in new_jobs:
    if (not job in previous_jobs):
      notifications.append("You submitted job " + job)
      sounds.append("Job submitted")

  for job in previous_jobs:
    if (not job in new_jobs):
      notifications.append("Job " + job + " ended")
      sounds.append("Job ended")
      continue	
    previous_status = previous_jobs[job]
    new_status = new_jobs[job]
    if (previous_status != new_status):
      if ("aswell" in new_status):
        notifications.append("Job " + job + "started")
        sounds.append("Job started")
  notify("\n".join(notifications))
  play_sounds(sounds)
  check_lt()

new_jobs = get_running_jobs(user_name, is_remote)
store_jobs(new_jobs)
while True:
	new_jobs = get_running_jobs(user_name, is_remote)
	previous_jobs = get_last_read_jobs()
	store_jobs(new_jobs)
	compare_jobs(previous_jobs, new_jobs)
	time.sleep(timer_interval)


