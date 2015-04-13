#!/bin/bash
######################################################
###### THIS SCRIPT IS MACHINE DEPENDENT so you'll have a different one in every repo
######################################################
# Now I'll do like this. this script must be copied in the main femus directory 
# and completed with the machine dependent stuff. then it must be run.

#problemi di questo script: se si interrompe non basta rigirarlo
# non basta che il file sia corrotto


#qual e' la cosa piu' importante? 
# E' fare in modo che TUTTI GLI ESEGUIBILI E LE LIBRERIE CHE SERVONO,
# in particolare MPI e PETSC, siano SEMPRE ACCESSIBILI,
# sia in fase di INSTALLAZIONE
# sia in fase di UTILIZZO!

# Vorrei fare questo in modo che MPI e PETSC siano utilizzabili
# SIA DA SOLI, SIA IN FUNZIONE DI ALTRE LIBRERIE come LA NOSTRA, LIBMESH, AND SO ON

# Ora, la questione e' che l'installazione di PETSC per FEMTTU e' DIVERSA 
# dall'installazione per FEMUS... pertanto, la configurazione per FEMTTU sara' diversa 
# dalla configurazione per FEMUS... 
# bisognerebbe fare un if sul CODICE anziche' sulla macchina...


################################################################################################
############## "Software" DIRECTORY ############
# Now we can install the software that we need #
################################################################################################
echo Install petsc

if test "$1" != "--prefix-external"; then
echo "The first argument must be --prefix-external"; exit;
fi

if test "$2" = ""; then
echo "The second argument must be the directory (either relative or absolute) where you want to install petsc"; exit;
fi


SOFTWARE_DIR=`readlink -f $2`
echo "=========" $SOFTWARE_DIR

mkdir -p $SOFTWARE_DIR
cd $SOFTWARE_DIR



#######################################################################
################# LA PROCEDURA UNIVERSALE ###########################
#############################################################


EXTRACT_COMMAND="tar xzvf"
TARGZ=".tar.gz"


#######################################################################

# We will take mpi directly from petsc

# echo Download, extract, compile an MPI implementation
# FM_MPI_DIR_REL=openmpi-1.6.5
# FM_MPI_DIR_ABS=$SOFTWARE_DIR/$FM_MPI_DIR_REL
# 
# # One should see if a certain installation already present is CORRECT... what is the definition of correct?!
# # Probably the petsc configure does something like that... for us it's too much, we just start over.
# echo =========== Remove previous installations
# rm -rf $FM_MPI_DIR_REL/
# if [ ! -f $FM_MPI_DIR_REL$TARGZ ]
# then
# echo =========== Download
# wget http://www.open-mpi.org/software/ompi/v1.6/downloads/$FM_MPI_DIR_REL$TARGZ   #where will he put this
# else
# echo The file $FM_MPI_DIR_REL$TARGZ already exists and we assume it is not corrupted
# fi
# echo =========== Extract
# $EXTRACT_COMMAND $FM_MPI_DIR_REL$TARGZ
# echo =========== Compile
# cd $FM_MPI_DIR_REL
# ./configure  --prefix=$PWD
# make -j all
# make install
# cd ..
# echo =========== Clean
# rm -f $FM_MPI_DIR_REL$TARGZ
# 
# 
# echo =========== Configure environment variables
# export PATH=$FM_MPI_DIR_ABS/bin:$PATH
# export LD_LIBRARY_PATH=$FM_MPI_DIR_ABS/lib64:$LD_LIBRARY_PATH
# 
# sleep 1


#######################################################################

# We will take hdf5 directly from petsc: he compiles the source. maybe it is not the last version but it's good that he compiles!

# # # # echo Download, extract, compile HDF5, or pre-binary, without compiling
# # # # FM_HDF5_DIR_REL=hdf5-1.8.12-linux-x86_64-shared
# # # # FM_HDF5_DIR_ABS=$SOFTWARE_DIR/$FM_HDF5_DIR_REL
# # # # 
# # # # 
# # # # echo =========== Remove previous installations
# # # # rm -rf $FM_HDF5_DIR_REL/
# # # # if [ ! -f $FM_HDF5_DIR_REL$TARGZ ]
# # # # then
# # # # echo =========== Download
# # # # wget http://www.hdfgroup.org/ftp/HDF5/current/bin/linux-x86_64/$FM_HDF5_DIR_REL$TARGZ
# # # # # wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.12.tar.gz
# # # # else
# # # # echo The file $FM_HDF5_DIR_REL$TARGZ already exists and we assume it is not corrupted
# # # # fi
# # # # echo =========== Extract
# # # # $EXTRACT_COMMAND $FM_HDF5_DIR_REL$TARGZ
# # # # echo =========== Clean
# # # # rm -f $FM_HDF5_DIR_REL$TARGZ
# # # # 
# # # # # source code
# # # # # wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.11.tar.gz
# # # # # now, before running configure, does it need to know the presence of mpi?
# # # # # here, i should put the flag to remove the valgrind error, and do maybe dbg and opt and else...
# # # # # ./configure --enable-using-memchecker --enable-clear-file-buffers
# # # # 
# # # # # wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.12.tar.gz
# # # # # extract 
# # # # # cd hdf5-1.8.12
# # # # # wget http://www.hdfgroup.org/ftp/lib-external/CMake/SZip.tar.gz
# # # # # wget http://www.hdfgroup.org/ftp/lib-external/CMake/ZLib.tar.gz
# # # # # wget http://www.hdfgroup.org/ftp/HDF5/examples/CMake/HDF518LinuxRWDICMake.cmake
# # # # # wget http://www.hdfgroup.org/ftp/HDF5/examples/CMake/CTestScript.cmake
# # # # # ctest -S HDF518LinuxRWDICMake.cmake,hdf5-1.8.12 -C RelWithDebInfo -VV -O hdf5.log
# # # # # 
# # # # # Now my question is: with this new system, how do we pass the CONFIGURE OPTIONS like we did before?
# # # # 
# # # # echo =========== Configure environment variables
# # # # export LD_LIBRARY_PATH=$FM_HDF5_DIR_ABS/lib:$LD_LIBRARY_PATH
# # # # 
# # # # # now, it seems like petsc does not find the libz which is automatically brought inside by hdf5... so let me set the library path for that
# # # # # the alternative is to download hdf5 but to put it OUTSIDE of petsc
# # # # # CHIARO CI VOGLIONO ANCHE GLI INCLUDE!!! E questi il pacchetto hdf5 NON CE LI HA!!! Quindi bisogna per forza prendere il pacchetto zlib da sistema o scaricarlo!!!
# # # # # Prendiamolo da sistema per ora
# # # # 
# # # # # HDF5 also needs sz... who is providing it? 
# # # # # Since HDF5 may need SZIP and ZLIB, i don't understand why they only put the PRECOMPILED LIBRARIES in the package WITHOUT the INCLUDE FILES...
# # # # # Well, if an application using HDF5 will not use those libraries DIRECTLY, then you don't need the includes but only the libs
# # # # # So... i have to do LD_LIBRARY_PATH to  let petsc reach libsz... now the only thing is that, doing this, i have TWO LIBZ, one from the system and one from hdf5.
# # # # # This may create a little conflict because if the SYSTEM LIBZ has the priority wrt the hdf5 package, then we link against a library which is not exactly
# # # # # the one with which hdf5 was compiled and linked... in this sense the most accurate solution would be to COMPILE HDF5... so far let us do with LD_LIBRARY_PATH
# # # # 
# # # # # The one thing i dont understand is why the petsc configure script DOES NOT FIND libz.so from the hdf5 lib folder. It should find it,
# # # # # because i added it in the LD_LIBRARY_PATH... WHY DOES IT NOT FIND IT?
# # # # 
# # # # # You need to set PATH and LD_LIBRARY_PATH  with the new MPI I think... at least PATH,
# # # # # otherwise where are the compilers?!?
# # # # # well, what you need to do is for sure to let petsc know about YOUR openmpi and YOUR HDF5...
# # # # # then i think that he will automatically find mpic++
# # # # # then, another point will be to let MY CODE be aware of MPI and HDF5,
# # # # # so I may need PATH and LD_LIBRARY_PATH for that
# # # # 
# # # # #plus, if what i download with petsc may be needed by future "higher level" packages i may consider installing them
# # # # # in a different folder than INSIDE petsc, for instance metis or parmetis, but we'll see next...
# # # # 
# # # # # Ok, the libraries are specified in LD_LIBRARY_PATH,
# # # # #     the binaries are specified in PATH,
# # # # # but how do I specify the INCLUDES?!?
# # # # # they must be specified in the COMPILING COMMAND IN THE MAKEFILE
# # # # # I think that if i installed hdf5 with petsc, i would not have any problem




sleep 1




#######################################################################
echo Download, extract, compile PETSC
FM_PETSC_DIR_REL=petsc
FM_PETSC_DIR_ABS=$SOFTWARE_DIR/$FM_PETSC_DIR_REL
export PETSC_DIR=$FM_PETSC_DIR_ABS

myarchs=(linux-opt linux-dbg)
debugflag=(0 1)


echo =========== Remove previous installations
rm -rf $FM_PETSC_DIR_REL/
echo =========== Clone
git clone -b maint https://bitbucket.org/petsc/petsc $FM_PETSC_DIR_REL

cd $FM_PETSC_DIR_REL

for i in 0 1
do
export PETSC_ARCH=${myarchs[i]}
echo =========== Configure
./configure  --with-debugging=${debugflag[i]}  --with-x=1  --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --with-clanguage=cxx  COPTFLAGS='-O3 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native' --with-shared-libraries=1 --download-openmpi=1 --download-fblaslapack=1 --download-blacs=1 --download-scalapack=1 --download-metis=1 --download-parmetis=1 --download-mumps=1 --download-hdf5=1
#   --with-mpi-dir=$FM_MPI_DIR_ABS 
#   --with-hdf5-dir=$FM_HDF5_DIR_ABS
echo =========== Make
make all #don't put -j here because it gives error! it's automatically -j!
make test
done
# I don't want to redownload everything once i have them... is it possible? i should tell it to download from the other linux-dbg folder... 
# or, when i download the external packages, i should put them OUTSIDE of the petsc folder...
# i think that the script is already prepared for that because the next time it was much faster

cd ..




sleep 1


#######################################################################
echo Download, extract, compile LIBMESH

# What we should do now.
# we should download libmesh with git
# then run the configure. Then tell to put libmesh in a different build directory
# in any case, some (a few of the files not in the build directory will be MODIFIED.
# First we make a commit. Then the next idea is to apply a PATCH to libmesh. 
# After I apply the patch I can start the configuration and so on.
# Now my only question is this: if i apply the patch to v.0.9, will i be able to apply it to v1.0, or future versions?
# That's a key thing. Let me see.

# Ora dovremo fare due installazioni di libmesh separate per le due versioni di petsc

# What to do to prepare a libmesh installation
# It's assumed that mpi and petsc (and what petsc needs... hdf5?!) are properly installed,
# So we only need to setup the environment for them

# Now, the question is: basically  libmesh needs for me mpi and petsc.
# Now petsc was installed with several optional libraries, among which hdf5.
# Now when i compile libmesh i get an error involving hdf5...
# Then I have to tell libmesh that hdf5 is there, i think...

# Now, the alternatives when you install libmesh in the configure script are:
# Tell libmesh directly with --with-thelibraryiwant-dir=,
# or modify PATH and   LD_LIBRARY_PATH...
# The best thing is always to use the options of the configure script...
# Plus, I don't want to use hdf5 with libmesh
FM_LIBMESH_DIR_REL=libmesh
FM_LIBMESH_BUILD=build
FM_LIBMESH_INSTALL=install

echo =========== Remove previous installations
rm   -rf $FM_LIBMESH_DIR_REL/
mkdir -p $FM_LIBMESH_DIR_REL/${FM_LIBMESH_BUILD}
echo =========== Download from git repository, my fork in github
git clone https://github.com/giorgiobornia/libmesh.git  ${FM_LIBMESH_DIR_REL}/${FM_LIBMESH_BUILD}
# git clone git://github.com/libMesh/libmesh.git ${FM_LIBMESH_DIR_REL}/${FM_LIBMESH_BUILD}
for i in 0 1
do
mkdir -p ${FM_LIBMESH_DIR_REL}/${FM_LIBMESH_INSTALL}-petsc-${myarchs[i]}
cd  ${FM_LIBMESH_DIR_REL}/${FM_LIBMESH_BUILD}
git checkout gambit_format
export PETSC_ARCH=${myarchs[i]}
echo =========== Compile
./configure --disable-cxx11 --with-methods="opt dbg" --prefix=$SOFTWARE_DIR/${FM_LIBMESH_DIR_REL}/${FM_LIBMESH_INSTALL}-petsc-${myarchs[i]}
# this script should be clever enough to find the external packages in petsc
make -j $(( $(nproc)-1 ))
make install
cd ../../
done



sleep 1
echo Download, extract, compile PRE and  POST UTILITIES:  paraview, hdfview
# how do i put in a command line something that is chosen with some alternatives with menus in .php page?!?
# wget http://www.hdfgroup.org/ftp/HDF5/hdf-java/hdfview/hdfview_install_linux64.bin
# sh hdfview_install_linux64.bin   #this will need graphical interface
sleep 1



############### TO INSTALL GAMBIT ##################
# You need to run the script, install openmotif22-libs (libXm.so.3)  # crappy alternative: ln -s libXm.so.4 libXm.so.3, WHICH DOES NOT WORK
# also, you have to install libstdc++33 (libstdc++.so.5)
# After that you have to put the license in the correct place
# Then it works
