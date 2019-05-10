install_path=/home/morelbt/github/BPP/install_bpp

mkdir ${install_path}
cd ${install_path}


git clone https://github.com/BioPP/bpp-core
git clone https://github.com/BioPP/bpp-seq
git clone https://github.com/BioPP/bpp-phyl

mkdir bpp-core-build
mkdir bpp-phyl-build
mkdir bpp-seq-build
mkdir install_bpp



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


