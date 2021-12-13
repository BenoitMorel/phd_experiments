from ete3 import SeqGroup

def read_msa(alignment_file):
  msa = None
  formats = ["fasta", "phylip_relaxed", "iphylip_relaxed", "phylip_interleaved"]
  for f in formats:
    try:
      msa = SeqGroup(alignment_file, f)
      return msa
    except:
      pass
  print("Cannot read msa " + alignment_file)
  return None 
