names="ofstream ifstream chrono vector cout cerr endl string min( map exception unordered_map istringstream runtime_error setprecision ostringstream shared_ptr fixed ostream"
files="core/*/* core/*/*/* JointSearch/* SpeciesRax/* GeneRax/*"
for file in $files
do
  echo $file
  sed -i "s/using namespace std\;//g" $file
  for name in $names
  do
    echo $name
    sed -i "s/$name\([^A-Za-z0-9]\)/std::$name\1/g" $file
  done
  sed -i "s/std::std::/std::/g" $file
  sed -i "s/\([A-Za-z0-9_\.][A-Za-z0-9_\.]*\)std::/\1/g" $file
  sed -i "s/unordered_map/std::unordered_map/g" $file
  sed -i "s/std::std::/std::/g" $file
  sed -i "s/include \([<\"]\)std::/include\ \1/g" $file
done 
