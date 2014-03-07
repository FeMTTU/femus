######### DIFFERENCE between WHAT THE LIBRARY HAS and WHAT EACH APPLICATION USES ##########
# the configuration of each application does not depend on the library
# every application decides for a small amount of things
# the fact of having or not MPI is a LIBRARY THING
# the fact of having or not PETSc is a LIBRARY THING
# and so on for HDF5 
# Also for libmesh we must change the philosophy: the library HAS LINKS WITH LIBMESH, then the fact
# that each application USES them or not is another story...
# Per ora, quello che ho io sono tante applicazioni, che usano TUTTE libfemus,
# e in aggiunta ce n'e' una che usa ANCHE LIBMESH
# Potrei farne una che usa ANCHE DEAL II, o FENICS, o cosi' via...

# L'idea e': fai HARD CODED in un FILE quello che e' fissato una volta per tutte (es, le librerie che HAI)
# Invece, fai IN MODO FLESSIBILE, cioe' con delle VARIABILI d'AMBIENTE tipicamente, quello che PUOI FACILMENTE CAMBIARE.




#######################
# If you put TWICE the same object in a linking command, then it says multiply defined!

#######OUT OF SOURCE BUILD #############
# La cosa importante se lanci l'eseguibile da un altro punto del file system 
# e' che l'eseguibile usi PERCORSI ASSOLUTI dei file!

######## HOW TO FIND A PACKAGE ###################
# cmake --help-module-list
# cmake --help-module FindBZip2. 
# Every cmake module introduces several variables that you can use

# Sembra che ci siano dei problemi nel trovare i pacchetti, come MPI o HDF5. 
# Per trovare HDF5 ho visto che lui cerca nei BINARI di HDF5, quindi devo aggiungere quelli nel PATH!!
# Ok, adesso gli ho fatto trovare i bin, sembra che non riesco a fargli trovar gli include...
# come fanno di solito gli script di configurazione a trovare la cartella di include? quello che e' di sistema e' PATH o LD_LIBRARY_PATH, 
# ma per il resto...

# I want to be able to locate the dirs of hdf5 without installing it as a system package but as a local package!!!

# se il findcmake di hdf5 non funziona, lo modifico io in modo che funzioni...

########### HARD-CODED CONFIGURE and RECOMPILING Software ################
# Se aggiorni il sistema operativo ad una nuova versione, credo che non ci sia bisogno di ricompilare 
# tutti i software che hai installato. A meno che non cambino gli header delle blas o di altre librerie base,
# dovrebbe andare tutto liscio cosi'.



 # A build type is just a bunch of flags for PREPROC, COMPILER, AND LINKER

 # - find the COMPILER (mpi)
# - depending on the MODE, add the appropriate COMPILER FLAGS
# - find all the EXTERNAL PACKAGES: MPI, then HDF5, then PETSC, then LIBMESH (these packages must have been INSTALLED IN ADVANCE!)
# - let all the make commands ECHO so i can see them
# I think i should do two projects, one for the library and one for the executable
# tell cmake to set the name as Make.common

# One project for the library
# One project for the executable



# COMPILING COMMAND
# mpicxx 
# -DNDEBUG -DFEMUS_HAVE_PETSC -DFEMUS_DETECTED_PETSC_VERSION_MAJOR=3 -DFEMUS_DETECTED_PETSC_VERSION_MINOR=4 -DFEMUS_DETECTED_PETSC_VERSION_SUBMINOR=0 -DHAVE_MPI -DPETSC_USE_DYNAMIC_LIBRARIES 
# -O2 -felide-constructors -fPIC -funroll-loops -fstrict-aliasing -Wdisabled-optimization  
# -I/home/gbornia/Software/femus/include/ 
# -I/home/gbornia/Software/femus/include/algebra/ 
# -I/home/gbornia/Software/femus/include/parallel/
# -I/home/gbornia/Software/femus/applications/temper/config/ 
# -I/home/gbornia/Software/femus/applications/temper/ 
# -I/home/gbornia/Software//hdf5-1.8.11-linux-x86_64-shared/include 
# -I/home/gbornia/Software//petsc-3.4.0/include 
# -I/home/gbornia/Software//petsc-3.4.0/bmake/linux-opt 
# -I/home/gbornia/Software//petsc-3.4.0/linux-opt/include     
# -c main.C 
# -o main.x86_64-unknown-linux-gnu.opt.o


# LINKING COMMAND
# mpicxx  
# -O2 -felide-constructors -fPIC -funroll-loops -fstrict-aliasing -Wdisabled-optimization
# -o temper-opt 
# /home/gbornia/Software/femus/applications/temper/lib/x86_64-unknown-linux-gnu_opt//libfemus.so 
# TempPhysics.x86_64-unknown-linux-gnu.opt.o TempQuantities.x86_64-unknown-linux-gnu.opt.o 
# TempBcIc.x86_64-unknown-linux-gnu.opt.o EqnNS.x86_64-unknown-linux-gnu.opt.o
# EqnT.x86_64-unknown-linux-gnu.opt.o main.x86_64-unknown-linux-gnu.opt.o  
# -L/home/gbornia/Software//hdf5-1.8.11-linux-x86_64-shared/lib 
# -lhdf5  -lsz -Wl,-rpath,/home/gbornia/Software//petsc-3.4.0/linux-opt/lib
# -L/home/gbornia/Software//petsc-3.4.0/linux-opt/lib  
# -lpetsc -Wl,-rpath,/home/gbornia/Software/petsc-3.4.0/linux-opt/lib -L/home/gbornia/Software/petsc-3.4.0/linux-opt/lib
# -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu_dist_3.3 -lsuperlu_4.3 
# -Wl,-rpath,/home/gbornia/Software//petsc-3.4.0/linux-opt/lib -lflapack -lfblas -lX11 -lparmetis -lmetis -lpthread 
# -Wl,-rpath,/home/gbornia/Software//hdf5-1.8.11-linux-x86_64-shared/lib 
# -L/home/gbornia/Software//hdf5-1.8.11-linux-x86_64-shared/lib -lhdf5_hl -lhdf5 -lz -Wl,-rpath,/home/gbornia/Software/openmpi-1.6.4/lib64
# -L/home/gbornia/Software/openmpi-1.6.4/lib64 -Wl,-rpath,/usr/lib64/gcc/x86_64-suse-linux/4.7 
# -L/usr/lib64/gcc/x86_64-suse-linux/4.7 -Wl,-rpath,/usr/x86_64-suse-linux/lib -L/usr/x86_64-suse-linux/lib 
# -lmpi_f90 -lmpi_f77 -lgfortran -lm -lgfortran -lm -lgfortran -lm -lgfortran -lm -lm -lm -lquadmath -lm -lmpi_cxx
# -lstdc++ -ldl -lmpi -lrt -lnsl -lutil -lgcc_s -lpthread -ldl 
# -Wl,-rpath,/home/gbornia/Software//hdf5-1.8.11-linux-x86_64-shared/lib -Wl,-rpath,/home/gbornia/Software//petsc-3.4.0/linux-opt/lib

# I have to tell him to compile the files in the src and in this application as well


# cmake, like autotools, is a BUILD SYSTEM
# The idea is: either i first try to compile some other cmake-based project,
# and then see how i can extend to my stuff
# The fact that the build may be OUT OF SOURCE is not too big a problem,
# only i have to pay attention WHERE the INPUT and OUTPUT files are located 
# with respect to the executable position
# the idea is to do a build directory with the examples SEPARATED, so build/ex1, build/ex2, et cetera.
# Then, where to put input and output? ... I guess it's better to put a build/ folder in every example


# Ok, I did a cmake for a very simple program
# The idea is to do this for a LIBRARY.
# One feature is: OUT OF SOURCE BUILDS
# Another feature is: it creates the makefile and the executable OUT OF SOURCE
# So far I do not want to do that. I want to create the library
# with a cmake-based common makefile
# Then, this makefile will be included by all the single applications
# my problem is that so far i CUSTOMIZE the LIBRARY for EACH APPLICATION.
# Then I'll do all the reasoning for a SINGLE APPLICATION, and then extend to arbitrary application

# posso passare ad un include file le opzioni da CMake, il quale ad esempio puo' leggere una variabile
# da 

# FIND COMPILER
# FIND OPENMPI
# FIND HDF5
# FIND PETSC
# FIND LIBMESH

# The idea is that at the end you will do make install which will have src, include, lib, bin and so on

# I have to do the cmake for the library, but here it is more for the executable.
# let us do it for the executable and then we move to the library only

# I have to create a BUILD directory for every application
# I want to find MPI as a compiler, not just gcc

#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8)

set(MYEXECUTABLE $ENV{FM_MYAPP}_$ENV{FM_FEMUS_METHOD})
set(FM_OBJ_SUFFIX _$ENV{FM_FEMUS_METHOD}.o)
set(FM_INCLUDE_DIR include)
set(FM_SOURCE_DIR src)
set(FM_ALGEBRA_DIR algebra)      #both for source and for include
set(FM_PARALLEL_DIR parallel)    #both for source and for include
set(FM_SOURCE_DIR src)
set(FM_APPS_DIR applications)
set(FM_APP_CONFIG_DIR config)
set(FM_APP_LIB_DIR lib)

set(MYLIBRARY femus_$ENV{FM_FEMUS_METHOD})


project(${MYEXECUTABLE})

############## EXECUTABLE ####################
#How to generate the list of SOURCE FILES with a WILDCARD
aux_source_directory(./ MYSRC)
add_executable(${MYEXECUTABLE}  ${MYSRC}) 
# i need to tell that all the object files have the suffix _method.o
# so far i could do three different build folders which contain all the objects related to that method
# when you create an executable, of course there is one .C file which is more important because it contains the main() function
# I think the linker detects it automatically...


################ LIBRARY #####################
aux_source_directory(../../${FM_SOURCE_DIR} MYLIBSRC0)
aux_source_directory(../../${FM_SOURCE_DIR}/${FM_ALGEBRA_DIR} MYLIBSRC1)
aux_source_directory(../../${FM_SOURCE_DIR}/${FM_PARALLEL_DIR} MYLIBSRC2)
add_library(${MYLIBRARY} SHARED ${MYLIBSRC0} ${MYLIBSRC1} ${MYLIBSRC2})   #TODO this is shared only if you are in opt or dbg, if you are prof it must be STATIC
#Where do you put the library?
#Where do you put the executable?

################ SHOULD I add_dependencies?! #####################


message("FEMuS is looking for MPI ================================")
# use _<lang>_, without it it is deprecated
find_package(MPI REQUIRED)
if(MPI_CXX_FOUND)
include_directories(${MPI_CXX_INCLUDE_PATH})
target_link_libraries(${MYEXECUTABLE} ${MPI_CXX_LIBRARIES})   # This one seems to duplicate the object files
target_link_libraries(${MYLIBRARY} ${MPI_CXX_LIBRARIES})   # This one seems to duplicate the object files
#  if(MPI_CXX_COMPILE_FLAGS) #perche' dovrei mettere questo if? Per velocizzare?
   set_target_properties(${MYEXECUTABLE} ${MYLIBRARY} PROPERTIES COMPILE_FLAGS "-DHAVE_MPI ${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DHAVE_MPI") #why does the other one not work!!!!!!!!!!
#  endif()
#  if(MPI_CXX_LINK_FLAGS)
   set_target_properties(${MYEXECUTABLE} ${MYLIBRARY} PROPERTIES LINK_FLAGS "-DHAVE_MPI ${MPI_CXX_LINK_FLAGS}")
#  endif()
endif(MPI_CXX_FOUND)


message("FEMuS is looking for HDF5 ================================")
find_package(HDF5 REQUIRED)   # QUIET, REQUIRED ,,..
if(HDF5_FOUND)
include_directories(${HDF5_INCLUDE_DIRS}) # NON E' HDF5_INCLUDE_PATH... devi guardare nella documentazione del modulo!!!
 target_link_libraries(${MYEXECUTABLE} ${HDF5_LIBRARIES})
 target_link_libraries(${MYLIBRARY} ${HDF5_LIBRARIES})
#  if(HDF5_COMPILE_FLAGS)
#    set_target_properties(${MYEXECUTABLE} ${MYLIBRARY} PROPERTIES COMPILE_FLAGS "${HDF5_COMPILE_FLAGS}")
# endif()
#  if(HDF5_LINK_FLAGS)
   set_target_properties(${MYEXECUTABLE} ${MYLIBRARY} PROPERTIES LINK_FLAGS "${HDF5_LINK_FLAGS}")
#  endif()
endif(HDF5_FOUND)


message("FEMuS is looking for PETSc ================================")

set(CMAKE_MODULE_PATH $ENV{FEMUS_DIR}/cmake/) #all my cmake modules are there...

find_package(PETSc)   #Can I put required as an argument here?!
if(PETSC_FOUND)
 include_directories(${PETSC_INCLUDES})
 target_link_libraries(${MYEXECUTABLE} ${PETSC_LIBRARIES})
 target_link_libraries(${MYLIBRARY} ${PETSC_LIBRARIES})
#  if(PETSC_DEFINITIONS)
   set_target_properties(${MYEXECUTABLE} ${MYLIBRARY} PROPERTIES COMPILE_FLAGS " ${PETSC_DEFINITIONS}  -DFEMUS_HAVE_PETSC -DFEMUS_DETECTED_PETSC_VERSION_MAJOR=${PETSC_VERSION_MAJOR} -DFEMUS_DETECTED_PETSC_VERSION_MINOR=${PETSC_VERSION_MINOR} -DFEMUS_DETECTED_PETSC_VERSION_SUBMINOR=${PETSC_VERSION_SUBMINOR}")
#  endif()
endif(PETSC_FOUND)
# Perche' se faccio message di queste variabili, la prima volta le stampa e le successive volte NO?!?
# Perche' certe variabili le fa vedere tutte le volte, e altre solo la prima?
# questo rerun di cmake mi sembra misterioso...
# so, IF I RERUN CMAKE MORE THAN ONCE THE VARIABLES DO NOT SEEM TO BE "PERSISTENT" ...
# But, if run once, it goes like i want!
message("${PETSC_INCLUDES}")   #tutte le volte
message("${PETSC_VERSION}")    #solo la prima
message("${PETSC_VERSION_MAJOR}") #solo la prima
message("${PETSC_VERSION_MINOR}") #solo la prima
message("${PETSC_VERSION_SUBMINOR}") #solo la prima



#
#   LIBMESH
#
# find_package(LIBMESH) #This is only for GENCASE


# include
# FEMuS
message("FEMuS is setting FEMuS ================================")
include_directories($ENV{FEMUS_DIR}/${FM_APPS_DIR}/$ENV{FM_MYAPP}/${FM_APP_CONFIG_DIR})
include_directories($ENV{FEMUS_DIR}/${FM_INCLUDE_DIR})
include_directories($ENV{FEMUS_DIR}/${FM_INCLUDE_DIR}/${FM_ALGEBRA_DIR})
include_directories($ENV{FEMUS_DIR}/${FM_INCLUDE_DIR}/${FM_PARALLEL_DIR})
# HDF5 NO MORE NEEDED # include_directories($ENV{FM_BASEPATH_TO_HDF5}/$ENV{FM_HDF5_FOLDER}/$ENV{FM_HDF5_INCLUDE})
# HDF5 NO MORE NEEDED # target_link_libraries(${MYEXECUTABLE} $ENV{FM_BASEPATH_TO_HDF5}/$ENV{FM_HDF5_FOLDER}/$ENV{FM_HDF5_LIB}/libhdf5.so)  #HDF5 puts these now
# HDF5 NO MORE NEEDED # target_link_libraries(${MYEXECUTABLE} $ENV{FM_BASEPATH_TO_HDF5}/$ENV{FM_HDF5_FOLDER}/$ENV{FM_HDF5_LIB}/libsz.so)    # HDF5 puts these now


#  ifeq ($(FM_FEMUS_METHOD),opt)
#   FEMUS_CPPFLAGS += -DNDEBUG
#   FEMUS_CXXFLAGS +=  -O2 -felide-constructors -fPIC -funroll-loops -fstrict-aliasing -Wdisabled-optimization
#   FEMUS_CFLAGS   += -O2 -fPIC -funroll-loops -fstrict-aliasing
# 
# endif
# ifeq ($(FM_FEMUS_METHOD),dbg)
#   FEMUS_CPPFLAGS += -DDEBUG  -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
#   FEMUS_CXXFLAGS +=  -O0 -felide-constructors -g -ansi -pedantic -W -Wall -Wno-long-long -Wunused -Wpointer-arith -Wformat -Wparentheses -fPIC -Woverloaded-virtual
#   FEMUS_CFLAGS   += -g  -Wimplicit -fPIC
# endif
# ifeq ($(FM_FEMUS_METHOD),pro)
#   FEMUS_CPPFLAGS += -DNDEBUG
#   FEMUS_CXXFLAGS +=  -O2 -felide-constructors -fPIC -funroll-loops -fstrict-aliasing -Wdisabled-optimization -pg
#   FEMUS_CFLAGS   += -O2 -fPIC -funroll-loops -fstrict-aliasing -pg 
#   FEMUS_LIBS += -pg
# # TODO also in the linking command?
#   fm-enable-shared = no
# endif


message("FEMuS is setting COMPILE FLAGS ================================")
 
if($ENV{FM_FEMUS_METHOD} STREQUAL opt)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DNDEBUG -O2 -felide-constructors -fPIC -funroll-loops -fstrict-aliasing -Wdisabled-optimization")
elseif($ENV{FM_FEMUS_METHOD} STREQUAL dbg)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DDEBUG -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC  -O0 -felide-constructors -g -ansi -pedantic -W -Wall -Wno-long-long -Wunused -Wpointer-arith -Wformat -Wparentheses -fPIC -Woverloaded-virtual" )
elseif($ENV{FM_FEMUS_METHOD} STREQUAL pro)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DNDEBUG -O2 -felide-constructors -fPIC -funroll-loops -fstrict-aliasing -Wdisabled-optimization -pg")
endif($ENV{FM_FEMUS_METHOD} STREQUAL opt)
  
# link
target_link_libraries(${MYEXECUTABLE} $ENV{FEMUS_DIR}/${FM_APPS_DIR}/$ENV{FM_MYAPP}/${FM_APP_LIB_DIR}/libfemus_$ENV{FM_FEMUS_METHOD}.so)
# The following correspond to all the source files BUT THE MAIN
# in fase di linking stanno aggiungendo DUE VOLTE TUTTA LA LISTA DELLE LIBRERIE, IL CHE NON MI PIACE PER NIENTE!
# Sembra che aggiungere due volte i .so non dia problemi, mentre aggiungere due volte i .o dia MULTIPLE DEFINITION!!!
# O sbaglio io a mettere i .o in target_link_libraries, ed e' possibile, pero' anche lui sbaglia a mettermi due volte la lista delle 
# librerie che non serve!!! anche perche' non so se e' la stessa cosa sia per i .so sia per i .a ...
# Qual e' la differenza tra add_library and target_link_libraries ?!? devo farne uno, o due, o entrambi?!?
# add_library e' come add_executable, E' PER IL TARGET!!!
foreach(loop_var EqnNS EqnT TempBcIc TempPhysics TempQuantities) 
#    target_link_libraries(${MYEXECUTABLE} $ENV{FEMUS_DIR}/${FM_APPS_DIR}/$ENV{FM_MYAPP}/${loop_var}${FM_OBJ_SUFFIX})
# add_library(${MYEXECUTABLE} $ENV{FEMUS_DIR}/${FM_APPS_DIR}/$ENV{FM_MYAPP}/${loop_var}.C)
endforeach(loop_var)

# IF YOU PUT AN OBJECT TWICE IT GIVES MULTIPLE DEFINITION!
#What if i want to specify a library with the path and the name separately?
# What if i want to add flags to the LINK command and not to the COMPILE command?
# how to put the -Wl-rpath things?

target_link_libraries(${MYEXECUTABLE} $ENV{FM_BASEPATH_TO_MPI}/$ENV{FM_MPI_FOLDER}/$ENV{FM_MPI_LIB}/libmpi_cxx.so)
target_link_libraries(${MYEXECUTABLE} $ENV{FM_BASEPATH_TO_MPI}/$ENV{FM_MPI_FOLDER}/$ENV{FM_MPI_LIB}/libmpi.so)
target_link_libraries(${MYEXECUTABLE} $ENV{FM_BASEPATH_TO_PETSC}/$ENV{FM_PETSC_FOLDER}/$ENV{PETSC_ARCH}/lib/libpetsc.so)



####### 
set(CMAKE_VERBOSE_MAKEFILE on)  #make VERBOSE=1
#perche' questo comando non funziona piu' dopo che introduco find_package MPI ?!? O anche altre volte, quando rimuovo tutto... forse lo devo mettere in fondo...
set(CMAKE_CXX_OUTPUT_EXTENSION  .o)  #perche' non funziona?!
set(CMAKE_C_OUTPUT_EXTENSION .o)     #perche' non funziona?!
#########################
#########################
#########################


# Ok now I have to tell him 
# to compile all the other sources 
# and to link to them and to libfemus.so

# add_library: i want to add it, but not in a separate folder... i just want to add separate .C files to link!

# I want to add 

# Now other open questions are: can i put the executable and all of it NOT IN A BUILD directory?
# Or, what happens when i do

# How to view environment variables?
# How to change the position of the executable?
# how to put compile flags?
# how to find the mpi compiler?
# do i have to use autoconf together with cmake?
# How to view the cmake variables?
# HOW DO I DECLARE THE DEPENDENCIES in the CMAKE FILE?
# I want to be able to perform all the operations that i perform in make with CMake as well!
# For cycles, patsubst and so on. I want to have the same things

# doing SEPARATE FOLDERS for SEPARATE BUILDS might be ok if i want to distinguish the BUILD MODES!
# Distinguish OBJECT FILES by ARCHITECTURE and LINUX VERSION, COMPILING MODE (HOSTTYPE)?!?

# The fact of doing out of source builds is maybe a way to distinguish
# DIFFERENT COMPILERS 
# DIFFERENT COMPILE MODES

# Whatever you want to modify in the COMPILATION (icc vs gcc, some flags other than others, ...) So, whatever compilation set you have,
# You can make SEPARATE BUILDS. They should be in SEPARATE DIRECTORIES: both the executable together with its common 

# For every compiler, I can build in SEVERAL MODES and that's it

# There is no point in doing make install for the single application.

# It doesnt matter too much where the executable is, as long as you specify the file paths with ABSOLUTE PATH!!

# What is the definition of ARCHITECTURE? 
# Can i include the compile mode in that? No, it is something more strictly connected to the HARDWARE.
# Now, since on a given computer we will only have one CPU, what we could switch in terms of "ARCHITECTURE" is the choice 
# of the compiler... or maybe, if i have a 64 bit CPU, i could switch between 32 mode and 64 mode...
# Even, the switch on the compiler might not be the case, in the sense that for a given CPU there should be ONE AD-HOC compiler...
# Well, actually the idea of gcc  is to have a HARDWARE INDEPENDENT compiler... so, one could switch between the gcc and a more ad-hoc compiler.
# What is important to be defined for an ARCHITECTURE? OPERATING SYSTEM, BITS, AND COMPILER?
# I need to understand what this "x86_64-unknown-linux-gnu" string is composed of

# Ok, so the idea is to concentrate the build files in ONE DIRECTORY:
# For us we have libfemus, which is composed of several objects on its own. Then, we have the .o in each application.

# I am afraid i should make SEVERAL CMAKELISTS, for the SOURCE FILES in the APPLICATION, and for the SOURCE FILES THAT FORM THE LIBRARY.
