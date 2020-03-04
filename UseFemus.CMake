# - UseFemus.CMake collects a set of CMake variables
# useful for building a vertical application outside of the Femus package

SET(FEMUS_LIBRARY_DIRS "/Users/anthonygruber/4Femus/MyFEMuS/lib64/")

#IF(UNIX)
SET(FEMUS_LIBRARIES "/usr/local/Cellar/petsc/3.12.4/lib/libpetsc.dylib;/usr/lib/liblapack.dylib;/usr/lib/libblas.dylib;/usr/local/Cellar/open-mpi/4.0.2/lib/libmpi_usempif08.dylib;/usr/local/Cellar/open-mpi/4.0.2/lib/libmpi_usempi_ignore_tkr.dylib;/usr/local/Cellar/open-mpi/4.0.2/lib/libmpi_mpifh.dylib;/usr/local/Cellar/open-mpi/4.0.2/lib/libmpi.dylib;/usr/local/Cellar/gcc/9.2.0_3/lib/gcc/9/libgfortran.dylib;/usr/local/Cellar/gcc/9.2.0_3/lib/gcc/9/libquadmath.dylib;/usr/lib/libm.dylib;/usr/lib/libc++.dylib;/usr/lib/libdl.dylib;/Users/anthonygruber/4Femus/MyFEMuS/lib64/libb64.a;/Users/anthonygruber/4Femus/MyFEMuS/lib64/libjsoncpp.dylib;/Users/anthonygruber/4Femus/MyFEMuS/lib64/libadept.a;FPARSER_LIBRARY-NOTFOUND;MPI_EXTRA_LIBRARY-NOTFOUND;HDF5_LIBRARIES-NOTFOUND;/Users/anthonygruber/4Femus/MyFEMuS/lib64/libfemus.${DYLIB};")
#ENDIF(UNIX)

#IF(APPLE)
#SET(FEMUS_LIBRARIES #"/usr/local/Cellar/petsc/3.12.4/lib/libpetsc.dylib;/usr/lib/liblapack.dylib;/usr/lib/libblas.dylib;/usr/local/Cellar/open-mpi/4.0.2/lib/libmpi_usempif08.dylib;/usr/local/Cellar/open-mpi/4.0.2/lib/libmpi_usempi_ignore_tkr.dylib;/usr/local/Cellar/open-mpi/4.0.2/lib/libmpi_mpifh.dylib;/usr/local/Cellar/open-mpi/4.0.2/lib/libmpi.dylib;/usr/local/Cellar/gcc/9.2.0_3/lib/gcc/9/libgfortran.dylib;/usr/local/Cellar/gcc/9.2.0_3/lib/gcc/9/libquadmath.dylib;/usr/lib/libm.dylib;/usr/lib/libc++.dylib;/usr/lib/libdl.dylib;/Users/anthonygruber/4Femus/MyFEMuS/lib64/libb64.a;/Users/anthonygruber/4Femus/MyFEMuS/lib64/libjsoncpp.dylib;/Users/anthonygruber/4Femus/MyFEMuS/lib64/libadept.a;FPARSER_LIBRARY-NOTFOUND;MPI_EXTRA_LIBRARY-NOTFOUND;HDF5_LIBRARIES-NOTFOUND;/Users/anthonygruber/4Femus/MyFEMuS/lib64/libfemus.${DYLIB};")
#ENDIF(APPLE)

SET(FEMUS_INCLUDES "/Users/anthonygruber/4Femus/MyFEMuS/external/adept/adept-1.1/include;/usr/local/include;/usr/local/include;FPARSER_INCLUDE_DIR-NOTFOUND;/Users/anthonygruber/4Femus/MyFEMuS/src/algebra;/Users/anthonygruber/4Femus/MyFEMuS/src/mesh;/Users/anthonygruber/4Femus/MyFEMuS/src/meshGencase;/Users/anthonygruber/4Femus/MyFEMuS/src/utils;/Users/anthonygruber/4Femus/MyFEMuS/src/quadrature;/Users/anthonygruber/4Femus/MyFEMuS/src/parallel;/Users/anthonygruber/4Femus/MyFEMuS/src/equations;/Users/anthonygruber/4Femus/MyFEMuS/src/solution;/Users/anthonygruber/4Femus/MyFEMuS/src/enums;/Users/anthonygruber/4Femus/MyFEMuS/src/fe;/Users/anthonygruber/4Femus/MyFEMuS/src/physics;/Users/anthonygruber/4Femus/MyFEMuS/include;")
