names="array make_shared ofstream ifstream chrono vector cout cerr endl string map exception unordered_map istringstream runtime_error setprecision ostringstream shared_ptr fixed ostream stack min( ios stringstream istreambuf_iterator numerical_limits unordered_set isnormal isnan numeric_limits hash"
files="core/* core/*/* core/*/*/* JointSearch/* SpeciesRax/* GeneRax/*"
for file in $files
do
  echo $file
  sed -i "s/using namespace std\;//g" $file
  for name in $names
  do
    sed -i "s/$name\([^A-Za-z0-9_]\)/std::$name\1/g" $file
    sed -i "s/$name$/std::$name/g" $file
  done
  sed -i "s/std::std::/std::/g" $file
  sed -i "s/\([A-Za-z0-9_\.][A-Za-z0-9_\.]*\)std::/\1/g" $file
  sed -i "s/unordered_map/std::unordered_map/g" $file
  sed -i "s/min(/std::min(/g" $file
  sed -i "s/max(/std::max(/g" $file
  sed -i "s/std::std::/std::/g" $file
  sed -i "s/include \([<\"]\)std::/include\ \1/g" $file
done 
