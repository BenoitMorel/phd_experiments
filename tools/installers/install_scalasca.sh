#!/bin/bash
install_home="$HOME/install"
install_prefix="$HOME/install/install_dir"
mkdir -p $install_home
cd $install_home
mkdir -p $install_prefix




# I had to run this on my personal computer
# I commented it out for centos (it's actually not required for phobos)
#apt install gfortran

install() {
  mkdir build
  cd build
  ../configure --prefix=$install_prefix
  make -j 10
  make install
  cd ..
}

install_scalasca() {
  wget http://apps.fz-juelich.de/scalasca/releases/scalasca/2.4/dist/scalasca-2.4.tar.gz
  tar -xzvf scalasca-2.4.tar.gz
  cd scalasca-2.4
  install
  cd ..
}

install_scorep() {
  wget https://www.vi-hps.org/cms/upload/packages/scorep/scorep-5.0-rc1.tar.gz
  tar -xzvf scorep-5.0-rc1.tar.gz
  cd scorep-5.0-rc1
  install
  cd ..
}

install_cube() {
  wget http://apps.fz-juelich.de/scalasca/releases/cube/4.4/dist/cubelib-4.4.2.tar.gz
  tar -xzvf cubelib-4.4.2.tar.gz
  cd cubelib-4.4.2
  install
  cd ..
}

install_cube_gui() {
  wget http://apps.fz-juelich.de/scalasca/releases/cube/4.4/dist/cubegui-4.4.2.tar.gz
  tar -xzvf cubegui-4.4.2.tar.gz
  cd cubegui-4.4.2
  # here, the configure command might need something like
  # ../configure -with-cubelib-config=$install_prefix/bin
  mkdir build
  cd build
  ../configure --prefix=$install_prefix -with-cubelib-config=$install_prefix/bin
  make -j 10
  make install
  cd ..
  cd ..
}

#install_scalasca
#install_scorep
#install_cube

if xhost >& /dev/null ; 
then 
  install_cube_gui
  echo ""
  echo ""
  echo "Found a X session. Scalasca GUI has been installed"
else
  echo ""
  echo ""
  echo "Did not find a X session. Scalasca GUI has NOT been installed"
  echo "If this machine does not support X, you can profile your programs on it, and then scp the output files to another machine, with the GUI installed. There is also a way to produce text profiling (see scalasca documentation) without the GUI"
fi

echo ""
echo "Please add the following directory to your path:"
echo "$install_prefix/bin"


