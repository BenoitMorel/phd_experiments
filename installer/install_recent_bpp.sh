cd ~/github/ALE_and_BPP/


git clone https://github.com/BioPP/bpp-core
git clone https://github.com/BioPP/bpp-seq
git clone https://github.com/BioPP/bpp-phyl

mkdir bpp-core-build
mkdir bpp-phyl-build
mkdir bpp-seq-build
mkdir install_bpp

install_path=/home/morelbt/github/BPP/install_bpp


cd bpp-core-build/
cmake ../bpp-core -DCMAKE_INSTALL_PREFIX=/usr/ -DBUILD_TESTING=FALSE -DCMAKE_INSTALL_PREFIX=${install_path}
make -j 40
make install

cd ../bpp-seq-build/
cmake ../bpp-seq -DCMAKE_INSTALL_PREFIX=/usr/ -DBUILD_TESTING=FALSE -DCMAKE_INSTALL_PREFIX=${install_path}
make -j 40
make install

cd ../bpp-phyl-build/
cmake ../bpp-phyl -DCMAKE_INSTALL_PREFIX=/usr/ -DBUILD_TESTING=FALSE -DCMAKE_INSTALL_PREFIX=${install_path}
make -j 40
make install

cd ..
cd ALE
mkdir build
cd build

cmake .. -DCMAKE_LIBRARY_PATH=${install_path}/lib64 -DCMAKE_INCLUDE_PATH=${install_path}/include


