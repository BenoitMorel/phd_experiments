import fast_rf_cells
import sys
import fam




def analyze(datadir, run):
  run_key = "true.true - " + run
  rf_cells = fast_rf_cells.load_rf_cells(datadir)
  per_family_score = []
  for family in rf_cells:
    score = rf_cells[family][run_key]
    rrf = float(score[0])/float(score[1])
    per_family_score.append((rrf, score, family))
  per_family_score.sort(key=lambda tup: tup[0]) 
  for tup in per_family_score:
    print(tup)


if __name__ == '__main__':
  if (len(sys.argv) < 2):
    print("Syntax: data_dir run")
    exit(1)
  print(" ".join(sys.argv))
  datadir = sys.argv[1]
  run = sys.argv[2]
  analyze(datadir, run)
