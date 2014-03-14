#!/bin/bash
######################################################
###### THIS SCRIPT IS MACHINE DEPENDENT so you'll have a different one in every repo
###### It will be run ONE LEVEL ABOVE the femus toplevel directory
######################################################

# This script must be used only once one knows that everything is properly installed...
# Se uno dei due vuole essere uno script che aggiorna permanentemente
# l'environment allora dev'essere quello chiamato ESTERNAMENTE
# Forse devo rassegnarmi al source, fare un mini-source solo 
# per mpi e petsc, e mettere tutto il resto HARD-CODED...
# Bisogna stare attenti perche' questo comunica con i makefile quindi 
# non si possono cambiare i nomi delle variabili senza cambiare anche nei makefile

# Variabili di ingresso a questo script:
# $FM_FEMUS_METHOD {opt dbg pro}  #don't need this anymore, everything is managed by CMAKE
# $FM_PETSC_METHOD {opt dbg}
# $FM_LIBMESH_METHOD {opt dbg}
########################

#maybe i should first run the configuration script to set a bunch of
# environment variables with SOURCE and then run the install thing,
# IN WHICH I SHOULD ALSO COMPILE FEMUS!!!

# Allora, io da un lato sto mettendo lo scaricamento e la compilazione.
# di tutte le librerie esterne
# Dall'altro, mi dovro' occupare della compilazione di FEMUS stesso
# QUASI TUTTE le cose potrebbero essere HARD-CODED in un qualche file per il makefile
# Pero' per MPI e PETSC ho bisogno di mettere qualcosa nelle environment variables...

# NEW PHILOSOPHY
# Within the SOFTWARE directory you put a script with the functions 
# to configure the most important packages so that they are usable


########

# I completely changed the structure of the script
#Now we source it to have the functions,
# then we call them separately

# Now with these return doesn't end the whole thing 
# but simply it goes to the next
# Dovrei sostituire i return con degli "abort",
# oppure fare dei controlli sul VALORE di RITORNO
# della funzione PRECEDENTE!


# hdf5:
# we might have the version from 
#     SYSTEM package
#     OUR OWN package
#     SALOME package...


# Also do with BOTH LIBMESH and FEMUS,
#  and maybe MPI, BLAS, LAPACK, HDF5... in PROFILE mode...
# At least, do we have "dbg" for the others?

#The negative thing when changing a SETUP (e.g. dbg-opt-opt or dbg-dbg-opt ...)
# is that you have to CLOSE THE SHELL and open another one...

# First install openmpi (you just need gcc, gfortran, g++,...)
# Then install blas and lapack
# Then, install Petsc (need MPI + blas + lapack, + other stuff also)
# Then, install Libmesh (either checkout from repository, or download the package...)
# First unpack and call it "libmesh...-petsc-opt"
# Do METHOD=opt and METHOD=dbg
#  Then unpack and call it "libmesh...-petsc-dbg"
# Do METHOD=opt and METHOD=dbg
# Then, install HDF5 (actually, it may be directly linked with Petsc also...)

#There are two main phylosophies:
# Use the system packages as much as possible
# Download the packages separately

# Problems: USER PACKAGES and HOME PACKAGES
# OpenMPI is already installed  in the system.
# In order to make it more important,
# always put it BEFORE in the PATH and LD_LIBRARY_PATH ...
# Actually, it is not so correct to REMOVE the 
# system installation because it is shared 
# by other programs, so it should remain there...
# you only have to assure that you give the PRIORITY
# to your own package rather than the system one.


# I should use this script for setting ONLY MPI, 
# or ONLY PETSC, so, the single libraries, if I want
# have to find a way to configure only MPI automatically
# if I don't say anything then pick all the libraries
# if I say "--only-mpi" then do only mpi


########## FEMuS configuration script ###########
# It is also good for setting the environment for 
# the installation of the external libraries.

#launch this script from where it is (because of $PWD)
#source configure_femus.sh

#which one?
#1) run "./configure_femus.sh" doesnt exit
#3) run "sh configure_femus.sh" doesnt exit
#2) run ". configure_femus.sh" exits
#4) run "source configure_femus.sh"  exits
# ...either 1 or 3
#1) doesnt export
#3) doesnt export
# those who do not exit do not even export,
# because they create a subshell i think.
# do not use exit in sourced scripts
#source the script, so that the exported
# variables remain in the current shell

# if test "$2" = ""; then 
#involve all the libraries 
   
# these are like GLOBAL variables after you do "source"
# 
export FM_EXTERNAL=../../external/
for FM_ARG in $*
do
 echo $FM_ARG
done
echo "==============================================================
Attention! Now you have to launch the function fm_set_femus from the top-level femus directory!!!
Therefore cd to that directory!!!! ==================";

##########################################
##########################################
##########################################
# This can be used also as independent MPI configuration
function fm_set_mpi() {

############# MACHINE DEPENDENT ###################
export FM_BASEPATH_TO_MPI=$PWD/$FM_EXTERNAL/
export FM_MPI_FOLDER=openmpi-1.6.5
export FM_MPI_BIN=bin
export FM_MPI_LIB=lib64
export FM_MPI_MAN=share/man
############# END MACHINE DEPENDENT ###################


###################
export FM_MPI_DIR=$FM_BASEPATH_TO_MPI/$FM_MPI_FOLDER
#these recursive commands are those for which you should not re-run this script
# you should just check if $FM_BASEPATH_TO_MPI/$FM_MPI_BIN is already in the PATH variable
# if [[ "$PATH" =~ "*$FM_BASEPATH_TO_MPI/$FM_MPI_BIN*" ]]
# then 
export PATH=$FM_MPI_DIR/$FM_MPI_BIN:$PATH  
# else
#    echo "Avoid rewriting"
# fi
   export LD_LIBRARY_PATH=$FM_MPI_DIR/$FM_MPI_LIB:$LD_LIBRARY_PATH
   # is this one really needed?!?
########################


######## MAN PAGES  ###########
export MANPATH=$FM_MPI_DIR/$FM_MPI_MAN:$MANPATH
################################


}


##########################################
##########################################
##########################################
# This function is only oriented to the FEMuS makefile
function fm_set_hdf5() {


############## MACHINE DEPENDENT ###################
#HDF5
#AAAA you may experience problems if you link against USER library 
# but you include SYSTEM headers, or vice versa
#the goal is to get both from the SAME part
#system package
# export FM_BASEPATH_TO_HDF5=/usr
# export FM_HDF5_FOLDER=
# export FM_HDF5_INCLUDE=include
# export FM_HDF5_LIB=lib64
#user package
export FM_BASEPATH_TO_HDF5=$PWD/$FM_EXTERNAL/
export FM_HDF5_FOLDER=hdf5-1.8.12-linux-x86_64-shared
#compiled version: hdf5-1.8.10-patch1/hdf5/
export FM_HDF5_INCLUDE=include
export FM_HDF5_LIB=lib
#compiled version: lib64
export FM_HDF5_BIN=bin
############## END MACHINE DEPENDENT ###################

export FM_HDF5_DIR=$FM_BASEPATH_TO_HDF5/$FM_HDF5_FOLDER
export PATH=$FM_HDF5_DIR/$FM_HDF5_BIN:$PATH  
# it seems like it doesnt find the includes, let us try adding LD_LIBRARY_PATH...
export LD_LIBRARY_PATH=$FM_HDF5_DIR/$FM_HDF5_LIB:$LD_LIBRARY_PATH  
export HDF5_ROOT=$FM_HDF5_DIR
# AAA questo e' fondamentale per poter far trovare hdf5!!! senza questo lo script FindHDF5.cmake non trova le HDF5_INCLUDE_DIRS
# if you dont set the PATH then he doesnt find HDF5_LIBRARIES HDF5_INCLUDE_DIRS; PATH is enough for finding the libraries,
# so LD_LIBRARY_PATH doesnt seem to be needed in this case
# If there are problems maybe it is the case to find some other FindHDF5.cmake on the web... then how would it be used?
# Se pero' hdf5 e' di sistema, allora il concetto di HDF_ROOT non esiste!
# HDF5_ROOT e' utile per cmake, PATH e LIBRARY PATH sono utili per TUTTO il SISTEMA
# Usare l'hint era pero' qualcosa che volevo evitare se possibile... in modo che il sistema trovasse hdf5 automaticamente...
# in realta' anche libmesh, per trovare i pacchetti ad esso esterni, a cui lui non contribuisce, ha bisogno di qualche hint.. 
# tipo PETSC: dir e arch... Quindi alla fine ci sta dai...
# Quello che magari non ci starebbe e' anche specificare i nomi delle cartelle include bin e lib... dovrebbe trovarsele da solo
# una volta che ha il root... questo e' un primo obiettivo dell'usare Cmake.
}



####################################################################################
####################################################################################
####################################################################################
# This function is oriented both to the FEMuS makefile and for PETSc standalone
# alright inside here you realize you have to 
function fm_set_petsc() {

#First you have to set mpi
fm_set_mpi

if test "$FM_PETSC_METHOD" = ""; then
  echo "Set the compiling mode for PETSc first: FM_PETSC_METHOD ";  return;
fi

############ MACHINE DEPENDENT ###################
export FM_BASEPATH_TO_PETSC=$PWD/$FM_EXTERNAL/
export FM_PETSC_FOLDER=petsc-3.4.3
export PETSC_ARCH=linux-$FM_PETSC_METHOD
############ END MACHINE DEPENDENT ###################

export PETSC_DIR=$FM_BASEPATH_TO_PETSC/$FM_PETSC_FOLDER


}





####################################################################################
####################################################################################
####################################################################################
# this is oriented to configuring LibMesh FOR Femus
function fm_set_libmesh() {

if test "$FM_PETSC_METHOD" = ""; then
  echo "Set the compiling mode for PETSc first: FM_PETSC_METHOD ";  return;
fi

if test "$FM_LIBMESH_METHOD" = ""; then
  echo "Set the compiling mode for LibMesh first: FM_LIBMESH_METHOD ";  return;
fi


############## MACHINE DEPENDENT ###################
FM_LIBMESH_DIR_REL=libmesh
FM_LIBMESH_INSTALL=install
export FM_BASEPATH_TO_LM=$PWD/$FM_EXTERNAL/
export FM_LM_FOLDER=$FM_LIBMESH_DIR_REL/$FM_LIBMESH_INSTALL-petsc-linux-$FM_PETSC_METHOD
############## END MACHINE DEPENDENT ###################


######## LIBMESH #######
########################
# conversion from the femus-libmesh wrapper variable to libmesh METHOD
  if test "$FM_LIBMESH_METHOD" = "opt"; then
export   METHOD=opt
elif test "$FM_LIBMESH_METHOD" = "dbg"; then
export   METHOD=dbg
elif test "$FM_LIBMESH_METHOD" = "pro"; then
export   METHOD=pro
fi

# === check ===
if test $FM_LIBMESH_METHOD != $METHOD; then
echo "Not correctly setting LibMesh compile mode"; return;
fi



}

##########################################
##########################################
##########################################


function fm_set_femus() {

#fm_set_mpi      # should be first because hdf5 may need it

# fm_set_hdf5     # let me remove this now, i'll take it from petsc so far

fm_set_petsc    #needs FM_PETSC_METHOD; fm_set_mpi is set within here

fm_set_libmesh  #needs FM_PETSC_METHOD and FM_LIBMESH_METHOD


######## FEMUS #########
export FEMUS_DIR=$PWD
########################




echo -e \
"============================================================
============================================================
================ Welcome to FEMuS =========================  
===== The method for FEMuS will be given by CMake \n
===== The method for LibMesh is $FM_LIBMESH_METHOD \n
===== The method for PETSc is  $FM_PETSC_METHOD \n
"

echo -e "======== BEWARE !!! ============
If the methods for FEMuS and for LibMesh differ,
you are very likely to encounter problems
 during compilation, due to different compile modes!
========================================="




}
