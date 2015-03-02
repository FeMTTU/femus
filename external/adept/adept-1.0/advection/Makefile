include ../Makefile_include
PROGRAM = run_advection_benchmark
CORE_OBJECTS = advection_schemes.o advection_schemes_AD.o
ADEPT_OBJECTS = advection_schemes_adept.o
CPPAD_OBJECTS = advection_schemes_cppad.o
ADOLC_OBJECTS = advection_schemes_adolc.o
SACADO_OBJECTS = advection_schemes_sacado.o
SACADO_FAD_OBJECTS = advection_schemes_sacado_fad.o

OBJECTS = $(CORE_OBJECTS) $(ADEPT_OBJECTS) $(CPPAD_OBJECTS) $(ADOLC_OBJECTS) $(SACADO_OBJECTS) $(SACADO_FAD_OBJECTS)

# If you don't have all the libraries installed then remove the
# unavailable ones here
LIBRARY_FLAGS = -DADEPT -DCPPAD -DADOLC -DSACADO -DSACADO_FAD

LIB_ADOLC = -ladolc
#LIB_ADOLC = -ladolc-tape-to-disk

# Edit to point at your own include directories
MY_INCLUDES = $(INCLUDES) -I$(HOME)/include -I/usr/local/include -I$(HOME)/src/trilinos-11.0.3

# Edit to point at your own library locations
MY_LIBS = $(LIBS)  -L$(HOME)/lib -L/usr/local/lib -L$(HOME)/src/trilinos-11.0.3 $(LIB_ADOLC) -lsacado -lteuchos

# DEFAULT TARGET (what happens when you type "make"):
# The program "run_advection_benchmark_reduced" uses only the Adept
# library, so doesn't require the others to be installed first
$(PROGRAM)_reduced: $(CORE_OBJECTS) $(ADEPT_OBJECTS) run_advection_benchmark.cpp
	$(CPP) $(CPPFLAGS) $(MY_INCLUDES) -DADEPT $(CORE_OBJECTS) $(ADEPT_OBJECTS) run_advection_benchmark.cpp $(MY_LIBS) -o $(PROGRAM)_reduced

# The full program "run_advection_benchmark" includes the capability
# to access all automatic differentiation libraries
$(PROGRAM): $(OBJECTS) run_advection_benchmark.o ../lib/libadept.a
	$(CPP) $(CPPFLAGS) $(OBJECTS) run_advection_benchmark.o $(MY_LIBS) -o $(PROGRAM)

# Type "make all" to get both versions of the program
all: $(PROGRAM)_reduced $(PROGRAM)

# The following programs use only Adept plus one of the other libraries
$(PROGRAM)_adolc: $(CORE_OBJECTS) $(ADOLC_OBJECTS) $(ADEPT_OBJECTS) run_advection_benchmark.cpp
	$(CPP) $(CPPFLAGS) $(MY_INCLUDES) -DADOLC -DADEPT $(CORE_OBJECTS) $(ADOLC_OBJECTS) $(ADEPT_OBJECTS) run_advection_benchmark.cpp $(MY_LIBS) -o $(PROGRAM)_adolc

$(PROGRAM)_cppad: $(CORE_OBJECTS) $(CPPAD_OBJECTS) $(ADEPT_OBJECTS) run_advection_benchmark.cpp
	$(CPP) $(CPPFLAGS) $(MY_INCLUDES) -DCPPAD -DADEPT $(CORE_OBJECTS) $(CPPAD_OBJECTS) $(ADEPT_OBJECTS) run_advection_benchmark.cpp $(MY_LIBS) -o $(PROGRAM)_cppad

$(PROGRAM)_sacado: $(CORE_OBJECTS) $(SACADO_OBJECTS) $(SACADO_FAD_OBJECTS) $(ADEPT_OBJECTS) run_advection_benchmark.cpp
	$(CPP) $(CPPFLAGS) $(MY_INCLUDES) -DSACADO -DSACADO_FAD -DADEPT $(CORE_OBJECTS) $(SACADO_OBJECTS) $(SACADO_FAD_OBJECTS) $(ADEPT_OBJECTS) run_advection_benchmark.cpp $(MY_LIBS) -o $(PROGRAM)_sacado

# The following are the rules for compiling different object files
run_advection_benchmark.o: run_advection_benchmark.cpp *.h ../include/adept.h
	$(CPP) $(CPPFLAGS) $(MY_INCLUDES) $(LIBRARY_FLAGS) -DNDEBUG -c $<

%.o: %.cpp *.h ../include/adept.h
	$(CPP) $(CPPFLAGS) $(MY_INCLUDES) -DNDEBUG -c $<

%_adept.o: %.cpp *.h ../include/adept.h
	$(CPP) $(CPPFLAGS) $(MY_INCLUDES) -DADEPT -c $< -o $@

%_cppad.o: %.cpp *.h
	$(CPP) $(CPPFLAGS) -DNDEBUG $(MY_INCLUDES) -DCPPAD -c $< -o $@

%_adolc.o: %.cpp *.h
	$(CPP) $(CPPFLAGS) $(MY_INCLUDES) -DADOLC -c $< -o $@

%_sacado.o: %.cpp *.h
	$(CPP) $(CPPFLAGS) $(MY_INCLUDES) -DSACADO -c $< -o $@

%_sacado_fad.o: %.cpp *.h
	$(CPP) $(CPPFLAGS) $(MY_INCLUDES) -DSACADO_FAD -c $< -o $@

clean:
	rm -f $(OBJECTS) $(PROGRAM) $(PROGRAM)_reduced run_advection_benchmark.o