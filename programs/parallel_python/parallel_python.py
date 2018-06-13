import os
import sys
import random
import math
import time
import concurrent.futures as futures
#import dispy
import mpi4py.futures as mpi
from mpi4py import MPI

def slow_function(iterations):
  res = 0.0
  for i in range(0, iterations):
    for j in range(0, 10000):
      res += math.exp(random.random())

  print(res)


def sequential_process(jobs):
  for job in jobs:
    slow_function(job)

def futures_process(jobs):
  executor = futures.ProcessPoolExecutor()
  for job in jobs:
    executor.submit(slow_function, job)
  executor.shutdown()


def dispy_process(jobs):
  print("disabled because of hang")
#  cluster = dispy.JobCluster(slow_function)
#  index = 0
#  dispy_jobs = []
#  for job in jobs:
#    dispy_job = cluster.submit(job) # it is sent to a node for executing 'compute'
#    dispy_job.id = index # store this object for later use
#    index += 1
#    dispy_jobs.append(dispy_job)
#  cluster.wait()
#  for dispy_job in dispy_jobs:
#    dispy_job()


def mpi_process(jobs, cores):
  with mpi.MPIPoolExecutor(cores) as executor:
    for job in jobs:
      executor.submit(slow_function, job)


def main_fct():

  if (len(sys.argv) != 5):
    print("syntax: python parallel_python.py job_size job_number implem cores")
    print("Implementations: sequential, futures, dispy, mpi")
    sys.exit(1)

  job_size = int(sys.argv[1])
  job_number = int(sys.argv[2])
  implem = sys.argv[3]
  cores = int(sys.argv[4])
  print("Job size: " + str(job_size))
  print("Job number: " + str(job_number))
  print("Implementation: " + implem)
  print("Cores: " + str(cores))

  start_time = time.time()
  jobs = [job_size] * job_number

  if ("sequential" == implem):
    sequential_process(jobs)
  elif ("futures" == implem):
    futures_process(jobs)
  elif ("mpi" == implem):
    mpi_process(jobs, cores)
  else:
    print("invalid implem")
    sys.exit(1)


  elapsed_time = time.time() - start_time
  print("elapsed time " + str(elapsed_time) + "s")

if __name__ == '__main__':
  main_fct()

