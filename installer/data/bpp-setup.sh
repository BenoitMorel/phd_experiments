#! /bin/bash
CORE_VERSION=2.2.0
SEQ_VERSION=2.2.0
PHYL_VERSION=2.2.0
POPGEN_VERSION=2.2.0
SEQOMICS_VERSION=2.2.0
PHYLOMICS_VERSION=2.2.0
RAA_VERSION=2.2.0
QT_VERSION=2.2.0
SUITE_VERSION=2.2.0
PHYVIEW_VERSION=0.4.0

PATH_GET=/usr/bin/wget
PATH_TAR=/usr/bin/tar
PATH_CMK=/hits/sw/shared/apps/CMake/3.5.2/bin/cmake
PATH_MAK=/usr/bin/make
PATH_INSTALL=${HOME}/install/bio++
BPP_DOWNLOAD_URL=http://biopp.univ-montp2.fr/repos/sources

clear
echo "+------------------------------------------------------------+"
echo "|        Bio++ Installation script, version 2.2.0            |"
echo "|   This script will download, compile and install Bio++.    |"
echo "+------------------------------------------------------------+"
echo ""

#--------------------------------------------------------------------
# Checking parameters:

#wget
until [ -f $PATH_GET ] && [ -x $PATH_GET ]
do
  echo "$PATH_GET does not exists or is not executable."
  echo -n "Full path toward wget command:"
  read PATH_GET
done
echo "wget command: $PATH_GET"
echo ""

#tar
until [ -f $PATH_TAR ] && [ -x $PATH_TAR ]
do
  echo "$PATH_TAR does not exists or is not executable."
  echo -n "Full path toward tar command:"
  read PATH_TAR
done
echo "tar command: $PATH_TAR"
echo ""

#cmake
until [ -f $PATH_CMK ] && [ -x $PATH_CMK ]
do
  echo "$PATH_CMK does not exists or is not executable."
  echo -n "Full path toward cmake command:"
  read PATH_CMK
done
echo "cmake command: $PATH_CMK"
echo ""

#make
until [ -f $PATH_MAK ] && [ -x $PATH_MAK ]
do
  echo "$PATH_MAK does not exists or is not executable."
  echo -n "Full path toward make command:"
  read PATH_MAK
done
echo "make command: $PATH_MAK"
echo ""

#install dir
echo -n "Installation directory [$PATH_INSTALL]:"
#read path
if [ "$path" != "" ] && [ "$path" != "y" ] && [ "$path" != "yes" ];
then
  until [ -d $path ]
  do 
    echo "Directory $path does not exist."
    echo -n "Installation directory:"
    read path
  done
  PATH_INSTALL=$path
fi

#Download all archives:

# Core:
CORE="yes"
echo -n "Do you want to download Bpp-Core [yes]/no?"
#read test
if [ "$test" != "yes" ] && [ "$test" != "y" ] && [ "$test" != "" ];
then
  echo "Ok, ending here!"
  exit
else
  $PATH_GET -nc $BPP_DOWNLOAD_URL/bpp-core-$CORE_VERSION.tar.gz
fi

if [ "$CORE" == "yes" ];
then
  # SeqLib:
  SEQ="yes"
  echo -n "Do you want to download Bpp-Seq [yes]/no?"
  read test
  if [ "$test" != "yes" ] && [ "$test" != "y" ] && [ "$test" != "" ];
  then
    SEQ="no"
  else
    $PATH_GET -nc $BPP_DOWNLOAD_URL/bpp-seq-$SEQ_VERSION.tar.gz
  fi
fi

if [ "$CORE" == "yes" ]  && [ "$SEQ" == "yes" ];
then
  # PhylLib:
  PHYL="yes"
  echo -n "Do you want to download Bpp-Phyl [yes]/no?"
  read test
  if [ "$test" != "yes" ] && [ "$test" != "y" ] && [ "$test" != "" ];
  then
    PHYL="no"
  else
    $PATH_GET -nc $BPP_DOWNLOAD_URL/bpp-phyl-$PHYL_VERSION.tar.gz
  fi
fi

if [ "$CORE" == "yes" ] && [ "$SEQ" == "yes" ];
then
  # PopGenLib:
  POPGEN="yes"
  echo -n "Do you want to download Bpp-PopGen [yes]/no?"
  read test
  if [ "$test" != "yes" ] && [ "$test" != "y" ] && [ "$test" != "" ];
  then
    POPGEN="no"
  else
    $PATH_GET -nc $BPP_DOWNLOAD_URL/bpp-popgen-$POPGEN_VERSION.tar.gz
  fi
fi

if [ "$CORE" == "yes" ] && [ "$SEQ" == "yes" ];
then
  # Bpp-Raa:
  RAA="yes"
  echo -n "Do you want to download Bpp-Raa [yes]/no?"
  read test
  if [ "$test" != "yes" ] && [ "$test" != "y" ] && [ "$test" != "" ];
  then
    RAA="no"
  else
    $PATH_GET -nc $BPP_DOWNLOAD_URL/bpp-raa-$RAA_VERSION.tar.gz
  fi
fi

if [ "$CORE" == "yes" ] && [ "$SEQ" == "yes" ] && [ "$PHYL" == "yes" ];
then
  # Bpp-Qt:
  QT="yes"
  echo -n "Do you want to download Bpp-Qt [yes]/no?"
  read test
  if [ "$test" != "yes" ] && [ "$test" != "y" ] && [ "$test" != "" ];
  then
    QT="no"
  else
    $PATH_GET -nc $BPP_DOWNLOAD_URL/bpp-qt-$QT_VERSION.tar.gz
  fi
fi

if [ "$CORE" == "yes" ] && [ "$SEQ" == "yes" ];
then
  # Bpp-Seq-Omics:
  SEQOMICS="yes"
  echo -n "Do you want to download Bpp-Seq-Omics [yes]/no?"
  read test
  if [ "$test" != "yes" ] && [ "$test" != "y" ] && [ "$test" != "" ];
  then
    SEQOMICS="no"
  else
    $PATH_GET -nc $BPP_DOWNLOAD_URL/bpp-seq-omics-$SEQOMICS_VERSION.tar.gz
  fi
fi

if [ "$CORE" == "yes" ] && [ "$SEQ" == "yes" ] && [ "$PHYL" == "yes" ];
then
  # Bpp-Phyl-Omics:
  PHYLOMICS="yes"
  echo -n "Do you want to download Bpp-Phyl-Omics [yes]/no?"
  read test
  if [ "$test" != "yes" ] && [ "$test" != "y" ] && [ "$test" != "" ];
  then
    PHYLOMICS="no"
  else
    $PATH_GET -nc $BPP_DOWNLOAD_URL/bpp-phyl-omics-$SEQOMICS_VERSION.tar.gz
  fi
fi

if [ "$CORE" == "yes" ] && [ "$SEQ" == "yes" ] && [ "$PHYL" == "yes" ];
then
  # BppSuite:
  SUITE="yes"
  echo -n "Do you want to download BppSuite [yes]/no?"
  read test
  if [ "$test" != "yes" ] && [ "$test" != "y" ] && [ "$test" != "" ];
  then
    SUITE="no"
  else
    $PATH_GET -nc $BPP_DOWNLOAD_URL/bppsuite/bppsuite-$SUITE_VERSION.tar.gz
  fi
fi

if [ "$CORE" == "yes" ] && [ "$SEQ" == "yes" ] && [ "$PHYL" == "yes" ];
then
  # BppPhyView:
  PHYVIEW="yes"
  echo -n "Do you want to download BppPhyView [yes]/no?"
  read test
  if [ "$test" != "yes" ] && [ "$test" != "y" ] && [ "$test" != "" ];
  then
    PHYVIEW="no"
  else
    $PATH_GET -nc $BPP_DOWNLOAD_URL/bppphyview/bppphyview-$PHYVIEW_VERSION.tar.gz
  fi
fi

# Now uncompress, configure, compile and install all libs

if [ "$CORE" == "yes" ];
then
  clear
  echo "+-----------------------------------------------------------+"
  echo "| Unpacking, configuring, building and installing Bpp-Core  |"
  echo "+-----------------------------------------------------------+"
  echo ""

  $PATH_TAR -xvzf bpp-core-$CORE_VERSION.tar.gz
  cd bpp-core-$CORE_VERSION
  $PATH_CMK -DCMAKE_INSTALL_PREFIX=$PATH_INSTALL .
  $PATH_MAK 
  $PATH_MAK install
  cd ..
fi

export CPPFLAGS=-I$PATH_INSTALL/include
export LDFLAGS=-L$PATH_INSTALL/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PATH_INSTALL/lib

if [ "$SEQ" == "yes" ];
then
  clear
  echo "+------------------------------------------------------------+"
  echo "| Unpacking, configuring, building and installing Bpp-Seq    |"
  echo "+------------------------------------------------------------+"
  echo ""
  $PATH_TAR -xvzf bpp-seq-$SEQ_VERSION.tar.gz
  cd bpp-seq-$SEQ_VERSION
  $PATH_CMK -DCMAKE_INSTALL_PREFIX=$PATH_INSTALL .
  $PATH_MAK 
  $PATH_MAK install
  cd ..
fi

if [ "$PHYL" == "yes" ];
then
  clear
  echo "+------------------------------------------------------------+"
  echo "| Unpacking, configuring, building and installing Bpp-Phyl   |"
  echo "+------------------------------------------------------------+"
  echo ""
  $PATH_TAR -xvzf bpp-phyl-$PHYL_VERSION.tar.gz
  cd bpp-phyl-$PHYL_VERSION
  $PATH_CMK -DCMAKE_INSTALL_PREFIX=$PATH_INSTALL .
  $PATH_MAK 
  $PATH_MAK install
  cd ..
fi

if [ "$POPGEN" == "yes" ];
then
  clear
  echo "+------------------------------------------------------------+"
  echo "| Unpacking, configuring, building and installing Bpp-PopGen |"
  echo "+------------------------------------------------------------+"
  echo ""
  $PATH_TAR -xvzf bpp-popgen-$POPGEN_VERSION.tar.gz
  cd bpp-popgen-$POPGEN_VERSION
  $PATH_CMK -DCMAKE_INSTALL_PREFIX=$PATH_INSTALL .
  $PATH_MAK 
  $PATH_MAK install
  cd ..
fi

if [ "$SEQOMICS" == "yes" ];
then
  clear
  echo "+---------------------------------------------------------------+"
  echo "| Unpacking, configuring, building and installing Bpp-Seq-Omics |"
  echo "+---------------------------------------------------------------+"
  echo ""
  $PATH_TAR -xvzf bpp-seq-omics-$SEQOMICS_VERSION.tar.gz
  cd bpp-seq-omics-$SEQOMICS_VERSION
  $PATH_CMK -DCMAKE_INSTALL_PREFIX=$PATH_INSTALL .
  $PATH_MAK 
  $PATH_MAK install
  cd ..
fi

if [ "$PHYLOMICS" == "yes" ];
then
  clear
  echo "+----------------------------------------------------------------+"
  echo "| Unpacking, configuring, building and installing Bpp-Phyl-Omics |"
  echo "+----------------------------------------------------------------+"
  echo ""
  $PATH_TAR -xvzf bpp-phyl-omics-$PHYLOMICS_VERSION.tar.gz
  cd bpp-phyl-omics-$PHYLOMICS_VERSION
  $PATH_CMK -DCMAKE_INSTALL_PREFIX=$PATH_INSTALL .
  $PATH_MAK 
  $PATH_MAK install
  cd ..
fi

if [ "$RAA" == "yes" ];
then
  clear
  echo "+------------------------------------------------------------+"
  echo "| Unpacking, configuring, building and installing Bpp-Raa    |"
  echo "+------------------------------------------------------------+"
  echo ""
  $PATH_TAR -xvzf bpp-raa-$RAA_VERSION.tar.gz
  cd bpp-raa-$RAA_VERSION
  $PATH_CMK -DCMAKE_INSTALL_PREFIX=$PATH_INSTALL .
  $PATH_MAK 
  $PATH_MAK install
  cd ..
fi

if [ "$QT" == "yes" ];
then
  clear
  echo "+------------------------------------------------------------+"
  echo "| Unpacking, configuring, building and installing Bpp-Qt     |"
  echo "+------------------------------------------------------------+"
  echo ""
  $PATH_TAR -xvzf bpp-qt-$QT_VERSION.tar.gz
  cd bpp-qt-$QT_VERSION
  $PATH_CMK -DCMAKE_INSTALL_PREFIX=$PATH_INSTALL .
  $PATH_MAK 
  $PATH_MAK install
  cd ..
fi

if [ "$SUITE" == "yes" ];
then
  clear
  echo "+------------------------------------------------------------+"
  echo "| Unpacking, configuring, building and installing BppSuite   |"
  echo "+------------------------------------------------------------+"
  echo ""
  $PATH_TAR -xvzf bppsuite-$SUITE_VERSION.tar.gz
  cd bppsuite-$SUITE_VERSION
  $PATH_CMK -DCMAKE_INSTALL_PREFIX=$PATH_INSTALL .
  $PATH_MAK 
  $PATH_MAK install
  cd ..
fi

if [ "$PHYVIEW" == "yes" ];
then
  clear
  echo "+------------------------------------------------------------+"
  echo "| Unpacking, configuring, building and installing BppPhyView |"
  echo "+------------------------------------------------------------+"
  echo ""
  $PATH_TAR -xvzf bppphyview-$PHYVIEW_VERSION.tar.gz
  cd bppphyview-$PHYVIEW_VERSION
  $PATH_CMK -DCMAKE_INSTALL_PREFIX=$PATH_INSTALL .
  $PATH_MAK 
  $PATH_MAK install
  cd ..
fi

