#!/bin/bash


echo -e "Source script for setting some environment variables for FEMuS to work\n"


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

 EXTERNAL_BASEPATH=`readlink -f $2`

############# MACHINE DEPENDENT ###################
 FM_BASEPATH_TO_MPI=$EXTERNAL_BASEPATH    #NOW ALL THESE VARIABLES will APPEAR IN THE ENVIRONMENT... even if you don't put EXPORT... TODO see if i can improve this
 FM_MPI_FOLDER=petsc/$PETSC_ARCH
 FM_MPI_BIN=bin
 FM_MPI_LIB=lib
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
   # is this one really needed?!? yes, to override system mpi installation
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

EXTERNAL_BASEPATH=`readlink -f $2`
PETSC_METHOD=$4


echo "The method for petsc is " $PETSC_METHOD

if [ "$PETSC_METHOD" != "opt" ] &&  [ "$PETSC_METHOD" != "dbg" ] ; then
  echo "Petsc method not supported";  return;
fi

############ MACHINE DEPENDENT ###################
       FM_BASEPATH_TO_PETSC=$EXTERNAL_BASEPATH
       FM_PETSC_FOLDER=petsc
export PETSC_ARCH=arch-linux2-cxx-$PETSC_METHOD
############ END MACHINE DEPENDENT ###################

export PETSC_DIR=$FM_BASEPATH_TO_PETSC/$FM_PETSC_FOLDER

# set mpi
fm_set_mpi --prefix-basepath $EXTERNAL_BASEPATH



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

EXTERNAL_BASEPATH=`readlink -f $2`
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
export FM_LM_FOLDER=$FM_LIBMESH_DIR_REL/$FM_LIBMESH_INSTALL-petsc-$PETSC_METHOD    #soon we'll avoid this export
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
echo " --prefix-external EXTERNAL_PATH --method-petsc {opt,dbg} --method-libmesh {opt,dbg,pro} ";
return;
fi

EXTERNAL_PREFIX=`readlink -f $2`
PETSC_METHOD=$4
LIBMESH_METHOD=$6

if test "$1" != "--prefix-external"; then
echo "Use --prefix-external for external packages"; return;
fi
if test "$3" != "--method-petsc"; then
echo "Use --method-petsc for Petsc"; return;
fi

# if test "$5" != "--method-libmesh"; then
# echo "Use --method-libmesh for Libmesh"; return;
# fi


#fm_set_mpi      # should be first because hdf5 may need it

# fm_set_hdf5     # let me remove this now, i'll take it from petsc so far

fm_set_petsc   --prefix-basepath $EXTERNAL_PREFIX --method  $PETSC_METHOD                 #needs FM_PETSC_METHOD; fm_set_mpi is set within here

# fm_set_libmesh --prefix-basepath $EXTERNAL_PREFIX --method  $LIBMESH_METHOD --method-petsc $PETSC_METHOD #needs FM_PETSC_METHOD and FM_LIBMESH_METHOD



echo -e \
"============================================================
============================================================
================ Welcome to FEMuS =========================  
===== The method for FEMuS will be given by CMake \n
===== The method for PETSc is  $PETSC_METHOD \n
"

# ===== The method for LibMesh is $LIBMESH_METHOD \n

# 
# echo -e "======== BEWARE !!! ============
# If the methods for FEMuS and for LibMesh differ,
# you are very likely to encounter problems
#  during compilation, due to different compile modes!
# ========================================="




}
