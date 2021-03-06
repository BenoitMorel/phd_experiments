import os
import sys
sys.path.insert(0, 'scripts')
import experiments as exp


def get_model(subst_model):
  return subst_model.split("+")[0]

def get_gamma_rates(subst_model):
  if ("G" in subst_model.split("+")):
    return 4
  else:
    return 1

def is_invariant_sites(subst_model):
  return "+I" in subst_model

def is_dna(subst_model):
  dna_models = ["JC", "GTR"]
  return get_model(subst_model) in dna_models

def get_phyldog_model(subst_model):
  return get_model(subst_model)

def get_raxml_model(subst_model):
  subst_model = subst_model.replace("+I", "+IC") # empirical invariant sites
  res = get_model(subst_model)
  if (res == "POISSON"):
    res = "PROTGTR{"
    res += "/".join(["1"] * 190)
    res += "}"
    res += "+FE"
  rest = subst_model.split("+")[1:]
  for r in rest:
    res += "+" + r
  return res

def get_mrbayes_preset_line(subst_model):
  if (get_model(subst_model) == "POISSON"):
    return "\tprset aamodelpr=fixed(poisson);\n"
  if (get_model(subst_model) == "LG"):
    return "\tprset aamodelpr=fixed(lg);\n"
  elif(get_model(subst_model) == "DAYHOFF"):
    return "\tprset aamodelpr=fixed(dayhoff);\n"
  elif(get_model(subst_model) == "WAG"):
    return "\tprset aamodelpr=fixed(wag);\n"
  elif(get_model(subst_model) == "JC"):
    return "\tprset statefreqpr=fixed(equal);\n"
  else:
    return ""

def get_deleterious_model(subst_model):
  if (get_model(subst_model) == "LG"):
    return "LG"
  elif(get_model(subst_model) == "WAG"):
    return "WAG"
  elif(get_model(subst_model) == "JTT"):
    return "JTT"
  elif(get_model(subst_model) == "POISSON"):
    return "UNIFORMAA"
  else:
    return "JC69"

def get_mrbayes_lset_line(subst_model):
  model = get_model(subst_model)
  line = "\t"
  rates = get_gamma_rates(subst_model)
  line += "lset "
  if (model == "GTR"):
    line += "nst=6"
  elif (model == "JC"):
    line += "nst=2"
  if (rates == 1):
    if (is_invariant_sites(subst_model)):
      line += " rates=propinv"
    else:
      line += " rates=equal"
  else:
    if (is_invariant_sites(subst_model)):
      line += " rates=invgamma ngammacat=" + str(rates)
    else:
      line += " rates=gamma ngammacat=" + str(rates)
  line += ";\n"
  return line

def extract_raxml_model(raxml_model_file):
  res = lambda:0
  line = open(raxml_model_file).readlines()[0]
  split = line.split("+")
  str1 = split[0]
  str2 = split[1]
  res.model = str1[0:str1.find("{")]
  res.rates = str1[str1.find("{")+1:str1.find("}")].split("/")
  res.frequencies = str2[str2.find("{")+1:str2.find("}")].split("/")
  return res

def build_default_dna_model():
  res = lambda:0
  res.model = "GTR"
  res.rates = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  res.frequencies = [0.25, 0.25, 0.25, 0.25]
  return res

def model_to_seqgen_cmd(model_obj):
  cmd = []
  cmd.append("-m")
  cmd.append(model_obj.model)
  cmd.append("-r")
  cmd.extend([str(i) for i in model_obj.rates])
  cmd.append("-f")
  cmd.extend([str(i) for i in model_obj.frequencies])
  return cmd



def get_model_samples(sample_name):
  samples = []
  if (sample_name == "dnadefault"):
    samples.append(build_default_dna_model())
  elif (sample_name == "dna4"):
    for sample_file in os.listdir(exp.dna4_model_samples):
      if (sample_file.endswith("bestModel")):
        sample_path = os.path.join(exp.dna4_model_samples, sample_file)
        samples.append(extract_raxml_model(sample_path))
  else:
    print("invalid sample name")
    return None
  return samples

def get_model_sample_names():
  return ["dnadefault", "dna4"]


if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax: python sequence_model.py model_file")
    exit(1)
  model = extract_model(sys.argv[1])
  print("model: " + str(model.model))
  print("rates: " + str(model.rates))
  print("frequencies: " + str(model.frequencies))

  seqgen_cmd = model_to_seqgen_cmd(model)
  print(" ".join(seqgen_cmd))

