# Include the Makefile defining compiler versions and options
include ../Makefile_include

# The -DADEPT_RECORDING_PAUSABLE option enables the pause_recording
# and continue_recording functionality and is used by test_adept,
# although it will run correctly (but slightly more slowly) without
# this flag
ADEPT_FLAGS = -DADEPT_RECORDING_PAUSABLE

# The objects to create
OBJECTS = algorithm.o algorithm_noad.o \
	test_adept.o test_adept_with_and_without_ad.o \
	test_radiances.o simulate_radiances.o
GSL_OBJECTS = test_gsl_interface.o state.o rosenbrock_banana_function.o

GSL_LIBS = -lgsl -lgslcblas

# Default target: compile first three programs (which only require
# Adept)
default: test_adept test_adept_with_and_without_ad test_radiances

# Type "make all" to also compile test_gsl_interface, which requires
# the GSL library
all: test_adept test_adept_with_and_without_ad test_radiances test_gsl_interface

# Test program 1
test_adept: $(OBJECTS) ../lib/libadept.a
	$(CPP) $(CPPFLAGS) $(ADEPT_FLAGS) algorithm.o test_adept.o $(LIBS) \
	-o test_adept

# Test program 2
test_adept_with_and_without_ad: $(OBJECTS) ../lib/libadept.a
	$(CPP) $(CPPFLAGS) $(ADEPT_FLAGS) algorithm.o algorithm_noad.o \
	test_adept_with_and_without_ad.o $(LIBS) \
	-o test_adept_with_and_without_ad

# Test program 3
test_radiances: $(OBJECTS) ../lib/libadept.a
	$(CPP) $(CPPFLAGS) $(ADEPT_FLAGS) simulate_radiances.o \
	test_radiances.o $(LIBS) -o test_radiances

# Test program 4
test_gsl_interface: $(GSL_OBJECTS) ../lib/libadept.a
	$(CPP) $(CPPFLAGS) $(ADEPT_FLAGS) $(GSL_OBJECTS) \
		$(GSL_LIBS) $(LIBS) -o test_gsl_interface

# The no-automatic-differentiation version of the algorithm: uses the
# -DADEPT_NO_AUTOMATIC_DIFFERENTIATION to produce a version of the
# algorithm that takes double rather than adouble arguments
algorithm_noad.o: algorithm.cpp *.h ../include/adept.h
	$(CPP) $(CPPFLAGS) $(ADEPT_FLAGS) $(INCLUDES) -c algorithm.cpp \
	-DADEPT_NO_AUTOMATIC_DIFFERENTIATION -o $@

# All other object files created by compiling the corresponding source
# file without this flag
%.o: %.cpp *.h ../include/adept.h
	$(CPP) $(CPPFLAGS) $(ADEPT_FLAGS) $(INCLUDES) -c $<

# Remove all object files and executables
clean:
	rm -f $(OBJECTS) $(GSL_OBJECTS) \
		test_adept test_adept_with_and_without_ad \
		test_radiances test_gsl_interface
