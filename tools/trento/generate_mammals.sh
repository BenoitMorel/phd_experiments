

initialsubsampling=0.05
species="CastorUUUCanadensis MusUUUCaroli MicrotusUUUOchrogaster MicrocebusUUUMurinus CebusUUUCapucinus AotusUUUNancymaae HomoUUUSapiens GorillaUUUGorilla MacacaUUUMulatta EquusUUUCaballus FelisUUUCatus CanisUUUFamiliaris LoxodontaUUUAfricana OrnithorhynchusUUUAnatinus"
mintaxa=15
maxtaxa=100

studentratio=0.1
studentreplicates=10


# get a subset of the families
python tools/families/generate_families_with_subsampling.py $bef/vertebrates188 $initialsubsampling 1

# only keep species of interest
dataset1=$bef/vertebrates188_subsample${initialsubsampling}_rep0
dataset2=$bef/vertebrates188_subsample${initialsubsampling}_rep0_pruned
python tools/families/generate_families_with_prunespecies.py $dataset1 $dataset2 raxml-ng GTR+G 1 0 $species

python tools/families/generate_families_with_multicopy_only.py $dataset2
dataset3=${dataset2}_multicopy
python tools/families/generate_families_with_taxon_number.py $dataset3 $mintaxa $maxtaxa
dataset4=${dataset3}_mintaxa${mintaxa}_maxtaxa${maxtaxa}

finaldataset=$bef/tuto_mammals
cp -rL $dataset4 $finaldataset

rm -rf $dataset1 $dataset2 $dataset3 $dataset4


# generate 10 student samples
python tools/families/generate_families_with_subsampling.py $bef/tuto_mammals $studentratio $studentreplicates

