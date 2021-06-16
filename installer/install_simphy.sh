#!/bin/bash

if [ -z "$1" ]
  then
    echo "No argument supplied (please specify user name)"
  fi
prefix=$1
temp=$prefix/temp
opt=$prefix/opt
rm -rf $temp
mkdir -p $prefix
mkdir -p $temp
mkdir -p $opt


cd $temp

wget http://www.sqlite.org/2015/sqlite-autoconf-3080900.tar.gz #SQlite3 sources
wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz #GSL sources
wget https://gmplib.org/download/gmp/gmp-6.0.0a.tar.bz2 #GMP sources
wget https://www.mpfr.org/mpfr-current/mpfr-4.1.0.tar.gz #MPFR sources
git clone https://github.com/adamallo/SimPhy --branch flexiblesim #SimPhy sources

tar -xvzf gsl*
tar -xvjf gmp*
tar -xvzf mpfr*
tar -xvzf sqlite*
rm *.gz *.bz2
mv gsl* gsl
mv gmp* gmp
mv mpfr* mpfr
mv sqlite* sqlite3


#step1: SQLite3 compilation and installation
cd $temp/sqlite3
./configure --prefix=$opt/local
make
make install
make clean
#Step2: GSL compilation and installation
cd $temp/gsl
./configure --prefix=$opt/local
make
make install
make clean
#Step3: GMP compilation and installation
cd $temp/gmp
./configure --prefix=$opt/local
make
make install
make clean

#Step4: Environment configuration

export C_INCLUDE_PATH=$opt/local/include
export LD_LIBRARY_PATH=$opt/local/lib
export LDFLAGS=-L$opt/local/lib
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'"$opt"'/local/lib' >> $HOME/.bashrc
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'"$opt"'/local/lib

#Step5: MPFR compilation and installation
cd $temp/mpfr
./configure --prefix=$opt/local --with-gmp=$opt/local
make
make install
make clean
#Step6: SimPhy compilation


cd $temp/SimPhy
make
make clean

rm -rf $prefix/SimPhy_1.0.2
mv $temp/SimPhy/ $prefix/SimPhy_1.0.2
echo "Binary in $prefix/SimPhy/bin"

