
datadir=$bef/tuto_mammals #_subsample0.1_rep8
model=GTR
cores=40

python tools/families/run_raxml_supportvalues.py $datadir $model 0 1 0 $cores 0


python scripts/generax/launch_speciesrax.py  $datadir random raxml-ng $model normald $cores
python tools/families/run_astral_pro.py $datadir/ raxml-ng $model
python tools/families/run_concatenation.py $datadir  min $model $cores
python tools/families/run_concatenation.py $datadir  max $model $cores

python tools/families/species_analyze.py $datadir


