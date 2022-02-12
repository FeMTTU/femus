#!/bin/bash

#you need to install Java Compiler Compiler, an OpenJDK development environment

# In opensuse, run as root (if you have another distro, figure out the equivalent needed packages)
# zypper in javacc
# zypper in java-1_7_0-openjdk-devel


#=====================
if test "$1" != "--prefix"; then
echo "The first argument must be --prefix"; exit;
fi

if test "$2" = ""; then
echo "The second argument must be the directory (either relative or absolute) where you want to install hdfjava"; exit;
fi


SOFTWARE_DIR=`readlink -f $2`
echo "=========" $SOFTWARE_DIR

mkdir -p $SOFTWARE_DIR
cd $SOFTWARE_DIR
#=====================



HDFJAVA_TAR=hdf-java-2.11.0.tar
HDFJAVA_CMAKE=HDFJAVALinuxCMake.cmake

mkdir -p hdfjava
cd hdfjava
wget http://www.hdfgroup.org/ftp/HDF5/hdf-java/current/cmake/SZip.tar.gz
wget http://www.hdfgroup.org/ftp/HDF5/hdf-java/current/cmake/ZLib.tar.gz
wget http://www.hdfgroup.org/ftp/HDF5/hdf-java/current/cmake/JPEG8b.tar.gz

wget http://www.hdfgroup.org/ftp/HDF5/hdf-java/current/cmake/HDF5.tar.gz
wget http://www.hdfgroup.org/ftp/HDF5/hdf-java/current/cmake/HDF4.tar.gz
wget http://www.hdfgroup.org/ftp/HDF5/hdf-java/current/cmake/CTestScript.cmake
wget http://www.hdfgroup.org/ftp/HDF5/hdf-java/current/cmake/${HDFJAVA_CMAKE}
wget http://www.hdfgroup.org/ftp/HDF5/hdf-java/current/src/${HDFJAVA_TAR}

ctest -S HDFJAVALinuxCMake.cmake -C Release -V -O hdf-java.log

cd hdf-java
cd build

echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
echo "!!!!!!!!!!! Answer the next two questions with  \"y\" and \"y\" !!!!!!!!!!!!!! "
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "

./HDFView-2.11.0-Linux.sh


sleep 5

cd HDFView-2.11.0-Linux/HDF_Group/HDFView/2.11.0

INSTALL_DIRECTORY=`pwd`

sed '/export INSTALLDIR/ c\ export INSTALLDIR='${INSTALL_DIRECTORY}' ' -i  bin/hdfview.sh
     
cp bin/hdfview.sh   $HOME/bin/hdfview
 