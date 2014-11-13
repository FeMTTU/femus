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

#launch this script from where it is (because of $PWD)  # Now i am REMOVING THIS CONSTRAINT, EVENTUALLY!!!
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


# # i=0
# # for FM_ARG in $*
# # do
# #  echo $FM_ARG
# # if (test $FM_ARG = "--prefix")
# # then
# # echo "Trovata" $i
# # fi
# # (( i += 1 ))
# # echo  ${i}
# # echo $1
# # done

echo  "Source script for setting the FEMuS environment"


while getopts ":a" opt; do
  case $opt in
    a)
      echo "-a was triggered, Parameter: $OPTARG" >&2
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
#       exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
#       exit 1
      ;;
  esac
done



##########################################
##########################################
##########################################
# This can be used also as independent MPI configuration
function fm_set_mpi() {
# $1 --prefix-basepath
# $2 value

if test "$1" = "--help"; then
echo " --prefix-basepath ABSOLUTE_PATH ";
return;
fi


if test "$1" != "--prefix-basepath"; then
echo "The first argument must be --prefix-basepath"; return 18;
fi

 EXTERNAL_BASEPATH=$2

############# MACHINE DEPENDENT ###################
 FM_BASEPATH_TO_MPI=$EXTERNAL_BASEPATH    #NOW ALL THESE VARIABLES will APPEAR IN THE ENVIRONMENT... even if you don't put EXPORT... TODO see if i can improve this
 FM_MPI_FOLDER=petsc/$PETSC_ARCH
 FM_MPI_BIN=bin
 FM_MPI_LIB=lib64
############# END MACHINE DEPENDENT ###################


###################
 FM_MPI_DIR=$FM_BASEPATH_TO_MPI/$FM_MPI_FOLDER
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

}



####################################################################################
####################################################################################
####################################################################################
# This function is oriented both to the FEMuS makefile and for PETSc standalone
# alright inside here you realize you have to 
function fm_set_petsc() {
# $1 --prefix-basepath
# $2 value
# $3 --method
# $4 value

if test "$1" = "--help"; then
echo " --prefix-basepath ABSOLUTE_PATH --method {opt,dbg} ";
return;
fi

if test "$1" != "--prefix-basepath"; then
echo "The first argument must be --prefix-basepath"; return 18;
fi

if test "$3" != "--method"; then
echo "The third argument must be --method"; return 18;
fi

EXTERNAL_BASEPATH=$2
PETSC_METHOD=$4


echo "The method for petsc is " $PETSC_METHOD


#First you have to set mpi
fm_set_mpi --prefix-basepath $EXTERNAL_BASEPATH


if [ "$PETSC_METHOD" != "opt" ] &&  [ "$PETSC_METHOD" != "dbg" ] ; then
  echo "Petsc method not supported";  return;
fi

############ MACHINE DEPENDENT ###################
       FM_BASEPATH_TO_PETSC=$EXTERNAL_BASEPATH
       FM_PETSC_FOLDER=petsc
export PETSC_ARCH=linux-$PETSC_METHOD
############ END MACHINE DEPENDENT ###################

export PETSC_DIR=$FM_BASEPATH_TO_PETSC/$FM_PETSC_FOLDER


}





####################################################################################
####################################################################################
####################################################################################
# this is oriented to configuring LibMesh FOR Femus
function fm_set_libmesh() {
# $1 --prefix-basepath
# $2 value
# $3 --method
# $4 value
# $5 --method-petsc
# $6 value

if test "$1" = "--help"; then
echo " --prefix-basepath ABSOLUTE_PATH --method {opt,dbg,pro} --method-petsc {opt,dbg}";
return;
fi




if test "$1" != "--prefix-basepath"; then
echo "Use --method "; return;
fi

if test "$3" != "--method"; then
echo "Use --method "; return;
fi

if test "$5" != "--method-petsc"; then
echo "Use --method-petsc for Petsc method "; return;
fi

EXTERNAL_BASEPATH=$2
LIBMESH_METHOD=$4
PETSC_METHOD=$6


if [ "$LIBMESH_METHOD" != "opt" ] &&  [ "$LIBMESH_METHOD" != "dbg" ] &&  [ "$LIBMESH_METHOD" != "pro" ] ; then
  echo "Libmesh method not supported";  return;
fi

if [ "$PETSC_METHOD" != "opt" ] &&  [ "$PETSC_METHOD" != "dbg" ] ; then
  echo "Petsc method not supported";  return;
fi



############## MACHINE DEPENDENT ###################
FM_LIBMESH_DIR_REL=libmesh
FM_LIBMESH_INSTALL=install
export FM_BASEPATH_TO_LM=$EXTERNAL_BASEPATH                                                       #soon we'll avoid this export 
export FM_LM_FOLDER=$FM_LIBMESH_DIR_REL/$FM_LIBMESH_INSTALL-petsc-linux-$PETSC_METHOD    #soon we'll avoid this export
############## END MACHINE DEPENDENT ###################


######## LIBMESH #######
########################
# conversion from the femus-libmesh wrapper variable to libmesh METHOD
  if test "$LIBMESH_METHOD" = "opt"; then
export   METHOD=opt
elif test "$LIBMESH_METHOD" = "dbg"; then
export   METHOD=dbg
elif test "$LIBMESH_METHOD" = "pro"; then
export   METHOD=pro
fi

# === check ===
if test $LIBMESH_METHOD != $METHOD; then
echo "Not correctly setting LibMesh compile mode"; return;
fi



}

##########################################
##########################################
##########################################


function fm_set_femus() {
#$1 --prefix                                             
#$2 full path of the prefix for femus                    
#$3 --prefix-external                                             
#$4 value
#$5 --method-petsc
#$6 values: {opt,dbg}
#$7 --method-libmesh
#$8  values: {opt,dbg}


if test "$1" = "--help"; then
echo " --prefix ABSOLUTE_PATH --prefix-external ABSOLUTE_PATH_2 --method-petsc {opt,dbg} --method-libmesh {opt,dbg,pro} ";
return;
fi

FEMUS_PREFIX=$2
EXTERNAL_PREFIX=$4
PETSC_METHOD=$6
LIBMESH_METHOD=$8

if test "$1" != "--prefix"; then
echo "Use --prefix for femus source dir"; return;
fi
if test "$3" != "--prefix-external"; then
echo "Use --prefix-external for external packages"; return;
fi
if test "$5" != "--method-petsc"; then
echo "Use --method-petsc for Petsc"; return;
fi
if test "$7" != "--method-libmesh"; then
echo "Use --method-libmesh for Libmesh"; return;
fi


#fm_set_mpi      # should be first because hdf5 may need it

# fm_set_hdf5     # let me remove this now, i'll take it from petsc so far

fm_set_petsc   --prefix-basepath $EXTERNAL_PREFIX --method  $PETSC_METHOD                 #needs FM_PETSC_METHOD; fm_set_mpi is set within here

fm_set_libmesh --prefix-basepath $EXTERNAL_PREFIX --method  $LIBMESH_METHOD --method-petsc $PETSC_METHOD #needs FM_PETSC_METHOD and FM_LIBMESH_METHOD


######## FEMUS #########
export FEMUS_DIR=$FEMUS_PREFIX
########################




echo -e \
"============================================================
============================================================
================ Welcome to FEMuS =========================  
===== The method for FEMuS will be given by CMake \n
===== The method for LibMesh is $LIBMESH_METHOD \n
===== The method for PETSc is  $PETSC_METHOD \n
"

echo -e "======== BEWARE !!! ============
If the methods for FEMuS and for LibMesh differ,
you are very likely to encounter problems
 during compilation, due to different compile modes!
========================================="




}
