import os

def get_elapsed_time_dico(dataset_dir):
  dico = {}
  try:
    with open(os.path.join(dataset_dir, "runtimes.txt")) as writer:
      for line in writer.readlines():
        split = line.split(" ")
        dico[split[0]] = split[1].replace("\n", "")
  except:
    pass
  return dico

def save_elapsed_time_dico(dataset_dir, dico):
  with open(os.path.join(dataset_dir, "runtimes.txt"), "w") as writer:
    for key in dico:
      writer.write(key + " " + dico[key] + "\n")

def save_elapsed_time(dataset_dir, method_key, elapsed_time):
  dico = get_elapsed_time_dico(dataset_dir)
  dico[method_key] = str(elapsed_time)
  save_elapsed_time_dico(dataset_dir, dico)




