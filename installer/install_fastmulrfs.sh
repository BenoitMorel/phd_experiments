cd ..

git clone https://github.com/ekmolloy/fastmulrfs.git
cd fastmulrfs/external

git clone https://github.com/pranjalv123/phylokit.git
cd phylokit/
mkdir build
cd build/
cmake -DCMAKE_INSTALL_PREFIX=../../install ../src
make install
cd ../../

git clone https://github.com/ekmolloy/phylonaut.git
cd phylonaut
mkdir build
cd build/
cmake -DCMAKE_INSTALL_PREFIX=../../install ../src
make install
cd ../../

git clone https://github.com/ekmolloy/FastRFS.git
cd FastRFS/
mkdir build
cd build/
cmake ../src
make

git clone https://github.com/smirarab/ASTRAL.git
cd ASTRAL
unzip Astral*zip
mv Astral ..



