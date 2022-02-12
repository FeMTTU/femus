# This Makefile creates an executable that can run the radar and lidar
# multiple scattering forward model described at
# http://www.met.reading.ac.uk/clouds/multiscatter/ but with the
# option to use different adjoint codes: the hand-coded ones (in
# *_AD.c) and those automatically produced by Adept, ADOL-C, CppAD and
# two versions of Sacado. This means that you need to install ADOL-C,
# CppAD and Sacado first. See the README file.

include ../Makefile_include

UTILITY_OBJECTS = config.o distribution.o utilities.o \
		dump_input.o range_spacing_is_regular.o print_algorithms.o
UNDIFFERENTIATED_OBJECTS = small_angle.o wide_angle_regrid.o
ADJOINT_OBJECTS = singlescatter_AD.o \
		wide_angle_AD.o tdts_AD.o fast_small_angle_AD.o \
		anisotropic_AD.o variance_AD.o
EXTRA_OBJECTS = fast_small_angle_lag.o jacobian.o jacobian_linear.o
CORE_OBJECTS = variance.o anisotropic.o \
		singlescatter.o multiscatter.o \
		wide_angle.o tdts.o fast_small_angle.o explicit.o simple_algorithm.o
ADEPT_OBJECTS = variance_adept.o anisotropic_adept.o \
		singlescatter_adept.o multiscatter_adept.o \
		wide_angle_adept.o tdts_adept.o fast_small_angle_adept.o \
		explicit_adept.o simple_algorithm_adept.o
CPPAD_OBJECTS = variance_cppad.o anisotropic_cppad.o \
		singlescatter_cppad.o multiscatter_cppad.o \
		wide_angle_cppad.o tdts_cppad.o fast_small_angle_cppad.o \
		explicit_cppad.o simple_algorithm_cppad.o
ADOLC_OBJECTS = variance_adolc.o anisotropic_adolc.o \
		singlescatter_adolc.o multiscatter_adolc.o \
		wide_angle_adolc.o tdts_adolc.o fast_small_angle_adolc.o \
		explicit_adolc.o simple_algorithm_adolc.o
SACADO_OBJECTS = variance_sacado.o anisotropic_sacado.o \
		singlescatter_sacado.o multiscatter_sacado.o \
		wide_angle_sacado.o tdts_sacado.o fast_small_angle_sacado.o \
		explicit_sacado.o simple_algorithm_sacado.o
SACADO_FAD_OBJECTS = variance_sacado_fad.o anisotropic_sacado_fad.o \
		singlescatter_sacado_fad.o multiscatter_sacado_fad.o \
		wide_angle_sacado_fad.o tdts_sacado_fad.o fast_small_angle_sacado_fad.o \
		explicit_sacado_fad.o simple_algorithm_sacado_fad.o

OBJECTS = $(UTILITY_OBJECTS) $(UNDIFFERENTIATED_OBJECTS) $(CORE_OBJECTS) \
	$(ADJOINT_OBJECTS) $(EXTRA_OBJECTS) $(ADEPT_OBJECTS) $(CPPAD_OBJECTS) \
	$(ADOLC_OBJECTS) $(SACADO_OBJECTS) $(SACADO_FAD_OBJECTS)
OBJECTS_REDUCED = $(UTILITY_OBJECTS) $(UNDIFFERENTIATED_OBJECTS) \
	$(CORE_OBJECTS) \
	$(ADJOINT_OBJECTS) $(EXTRA_OBJECTS) $(ADEPT_OBJECTS)

LIBMS = ./libmultiscatter.a
LIBMS_REDUCED = ./libmultiscatter_reduced.a

PROGRAM = ./multiscatter
PROGRAM_REDUCED = ./multiscatter_reduced

#CPPFLAGS = -O3

# You will need to point the following to the locations of the ADOL-C,
# CppAD and Sacado header files:
MS_INCLUDES = $(INCLUDES) -I$(HOME)/include -I/usr/local/include -I$(HOME)/src/trilinos-11.0.3

# You will need to point the following to the location of the ADOL-C
# shared library:
MS_LIBS_REDUCED = $(LIBS) -L.
MS_LIBS = $(MS_LIBS_REDUCED) -L$(HOME)/lib -L/usr/local/lib -L$(HOME)/src/trilinos-11.0.3 -ladolc -lsacado -lteuchos


# Default target: compile only "multiscatter_reduced", which only uses
# Adept
default: $(PROGRAM_REDUCED)

# "make all" will compile the multiscatter library and then the full
# multiscatter program that uses all the automatic differentiation
# libraries
all: $(PROGRAM_REDUCED) $(PROGRAM)

# Full multiscatter program including all the automatic
# differentiation libraries
$(PROGRAM): $(LIBMS) multiscatter_ascii.o ../lib/libadept.a
	$(CPP) $(CPPFLAGS) $(MS_INCLUDES) multiscatter_ascii.o \
		-L../lib -lmultiscatter $(MS_LIBS) -o $(PROGRAM)

# The "reduced" version of the program includes only the original
# algorithm, the hand-coded and Adept versions. This is because on
# some platforms, having all the other libraries compiled into one
# executable slows down Adept.
$(PROGRAM_REDUCED): $(LIBMS_REDUCED) multiscatter_ascii.o ../lib/libadept.a
	$(CPP) $(CPPFLAGS) $(MS_INCLUDES) multiscatter_ascii.o \
		-L../lib -lmultiscatter_reduced $(MS_LIBS_REDUCED) \
		-o $(PROGRAM_REDUCED)

$(LIBMS): $(OBJECTS) multiscatter_AD.o
	ar -r $(LIBMS) $(OBJECTS) multiscatter_AD.o

$(LIBMS_REDUCED): $(OBJECTS_REDUCED) multiscatter_reduced_AD.o
	ar -r $(LIBMS_REDUCED) $(OBJECTS_REDUCED) multiscatter_reduced_AD.o

# Each of the main source files is compiled six times, with different
# -Dxxx compiler directives to generate different object files. With
# no extra compiler directive the original forward algorithm is
# produced. With -DADEPT, -DCPPAD, -DADOLC, etc, the different automatic
# differentiation tools are applied (see ms_switch.h for how). They
# produce functions overloaded for the different active types of these
# three libraries. The -DNDEBUG flag doubles the speed of CppAD but
# omitting error checking.
%.o: %.c ms.h ms_switch.h multiscatter.h
	$(CPP) $(CPPFLAGS) -DNDEBUG $(MS_INCLUDES) -c $<

%_adept.o: %.c ms.h ms_switch.h multiscatter.h ../include/adept.h
	$(CPP) $(CPPFLAGS) $(MS_INCLUDES) -DADEPT -c $< -o $@

%_cppad.o: %.c ms.h ms_switch.h multiscatter.h
	$(CPP) $(CPPFLAGS) -DNDEBUG $(MS_INCLUDES) -DCPPAD -c $< -o $@

%_adolc.o: %.c ms.h ms_switch.h multiscatter.h
	$(CPP) $(CPPFLAGS) $(MS_INCLUDES) -DADOLC -c $< -o $@

%_sacado.o: %.c ms.h ms_switch.h multiscatter.h
	$(CPP) $(CPPFLAGS) $(MS_INCLUDES) -DSACADO -c $< -o $@

%_sacado_fad.o: %.c ms.h ms_switch.h multiscatter.h
	$(CPP) $(CPPFLAGS) $(MS_INCLUDES) -DSACADO_FAD -c $< -o $@

multiscatter_AD.o: ms.h ms_switch.h multiscatter.h ../include/adept.h
	$(CPP) $(CPPFLAGS) -DNDEBUG $(MS_INCLUDES) -c multiscatter_AD.c

multiscatter_reduced_AD.o: ms.h ms_switch.h multiscatter.h ../include/adept.h
	$(CPP) $(CPPFLAGS) -DSKIP_NON_ADEPT -DNDEBUG $(MS_INCLUDES) \
		-c multiscatter_AD.c -o multiscatter_reduced_AD.o

clean:
	rm -f $(OBJECTS) $(LIBMS) $(LIBMS_REDUCED) \
		$(PROGRAM) $(PROGRAM_REDUCED) multiscatter_ascii.o \
		multiscatter_AD.o multiscatter_reduced_AD.o

clean-autosaves:
	rm -f *~
