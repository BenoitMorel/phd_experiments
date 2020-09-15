
if [ -z "$1" ]
  then
    echo "No argument supplied (please specify user name)"
  fi
user=$1

p=/home/$user/temp

#: '

mkdir -p $p
mkdir -p /home/$user/opt
cd $p
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
cd $p/sqlite3
./configure --prefix=/home/$user/opt/local
make
make install
make clean
#Step2: GSL compilation and installation
cd $p/gsl
./configure --prefix=/home/$user/opt/local
make
make install
make clean
#Step3: GMP compilation and installation
cd $p/gmp
./configure --prefix=/home/$user/opt/local
make
make install
make clean

#Step4: Environment configuration

export C_INCLUDE_PATH=/home/$user/opt/local/include
export LD_LIBRARY_PATH=/home/$user/opt/local/lib
export LDFLAGS=-L/home/$user/opt/local/lib
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/'"$user"'/opt/local/lib' >> /home/$user/.bashrc

#Step5: MPFR compilation and installation
cd /home/$user/temp/mpfr
./configure --prefix=/home/$user/opt/local --with-gmp=/home/$user/opt/local
make
make install
make clean
#Step6: SimPhy compilation


cd /home/$user/temp/SimPhy
make
make clean

echo "Binary in $p/SimPhy/bin"


