import os

def get_metrics(dataset_dir, metric_name):
  dico = {}
  try:
    with open(os.path.join(dataset_dir, "metrics", metric_name + ".txt")) as writer:
      for line in writer.readlines():
        split = line.split(" ")
        dico[split[0].lower()] = split[1].replace("\n", "")
  except:
    return None
  return dico

def get_metrics_methods(dataset_dir, metric_name):
  methods = []
  try:
    with open(os.path.join(dataset_dir, "metrics", metric_name + ".txt")) as writer:
      for line in writer.readlines():
        split = line.split(" ")
        methods.append(split[0])
  except:
    pass
  return methods

def save_dico(dataset_dir, dico, metric_name):
  try:
    os.makedirs(os.path.join(dataset_dir, "metrics"))
  except:
    pass
  with open(os.path.join(dataset_dir, "metrics", metric_name + ".txt"), "w") as writer:
    for key, value in sorted(dico.items(), key=lambda x: float(x[1])):
      writer.write(key + " " + value + "\n")

def save_metrics(dataset_dir, method_key, metric_value, metric_name):
  method_key = method_key.lower()
  dico = get_metrics(dataset_dir,  metric_name)
  if (dico == None):
    dico = {}
  dico[method_key] = str(metric_value)
  save_dico(dataset_dir, dico, metric_name)




