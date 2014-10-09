# - Find the Function Parser for C++ library
# This C++ library offers a class which can be used to parse and evaluate a mathematical function 
# from a string (which might be eg. requested from the user). The syntax of the function string is 
# similar to mathematical expressions written in C/C++.
# The function can then be evaluated with different values of variables.
# For example, a function like "sin(sqrt(x*x+y*y))" can be parsed from a string 
# (either std::string or a C-style string) and then evaluated with different values of x and y. 
# This library can be useful for evaluating user-inputted functions, or in some cases 
# interpreting mathematical expressions in a scripting language.
# This library aims for maximum speed in both parsing and evaluation, while keeping maximum portability. 
# The library should compile and work with any standard-conforming C++ compiler.
# Different numerical types are supported: double, float, long double, long int, 
# std::complex (of types double, float and long double), multiple-precision floating point 
# numbers using the MPFR library, and arbitrary precision integers using the GMP library. 
# (Note that it's not necessary for these two libraries to exist in the system in order 
# to use the Function Parser library with the other numerical types. 
# Support for these libraries is optionally compiled in using preprocessor settings.)
# This Library is distributed under the Lesser General Public License (LGPL) version 3. 
#
# === Variables ===
#
# FPARSER_FOUND               - True if found, otherwise all other vars are undefined
# FPARSER_INCLUDE_DIR         - The include dir for main *.h files
# FPARSER_LIBRARY             - The library to link against
#
#

INCLUDE (SelectLibraryConfigurations)
INCLUDE (FindPackageHandleStandardArgs)

SET(FPARSER_LIBRARY FPARSER_LIBRARY-NOTFOUND CACHE STRING "FParser library to link against" FORCE)
SET(FPARSER_INCLUDE_DIR FPARSER_INCLUDE_DIR-NOTFOUND CACHE STRING "Fparser include directory" FORCE )

# The HINTS option should only be used for values computed from the system.
SET (_FPARSER_HINTS
    $ENV{HOME}/.local
)
# Hard-coded guesses should still go in PATHS. This ensures that the user
# environment can always override hard guesses.
SET (_FPARSER_PATHS
    $ENV{HOME}/.local
    /usr/lib64
    /usr
    /usr/include
)

FIND_LIBRARY (FPARSER_LIBRARY "libfparser.so"
    HINTS ${_FPARSER_HINTS}
    PATHS ${_FPARSER_PATHS}
    PATH_SUFFIXES 
)

FIND_PATH (FPARSER_INCLUDE_DIR "fparser.hh"
    HINTS ${_FPARSER_HINTS}
    PATHS ${_FPARSER_PATHS}
    PATH_SUFFIXES 
)

SET ( FPARSER_INCLUDE_DIRS ${FPARSER_INCLUDE_DIR})
SET ( FPARSER_LIBRARIES ${FPARSER_LIBRARY})

SET (FPARSER_FOUND "FALSE")
IF (FPARSER_LIBRARY)
  SET (FPARSER_FOUND "TRUE")
ENDIF (FPARSER_LIBRARY)

