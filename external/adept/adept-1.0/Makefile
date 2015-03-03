# Type "make" to make the Adept library and put it in the "lib"
# directory, then compile the test files in "test". This will also
# compile the Adept-only benchmark programs
# advection/run_advection_benchmark_reduced and
# multiscatter/multiscatter_reduced.

# If you have the ADOL-C, Sacado and CppAD libraries installed on your
# system, type "make all" to try to compile the programs used to
# produce the benchmarks in the Adept paper:
# advection/run_advection_benchmark and multiscatter/multiscatter.

# You also have the "make install" and "make clean" which do what
# you'd expect.

# Note that the Makefiles in subdirectories include the
# Makefile_include file in this directory, which includes details
# about compilers and compiler flags

PREFIX = /usr/local
LIBDIR = $(PREFIX)/lib
INCLUDEDIR = $(PREFIX)/include

default: libadept.a test advection multiscatter

all: default advection-all multiscatter-all

libadept.a:
	cd adept && $(MAKE)

test: libadept.a
	cd test && $(MAKE)

advection: libadept.a
	cd advection && $(MAKE)

advection-all: libadept.a
	cd advection && $(MAKE) all

multiscatter:
	cd multiscatter && $(MAKE)

multiscatter-all:
	cd multiscatter && $(MAKE) all

install: libadept.a
	mkdir -p $(INCLUDEDIR)
	install -m 0644 include/adept.h $(INCLUDEDIR)
	mkdir -p $(LIBDIR)
	install -m 0644 lib/libadept.a $(LIBDIR)

uninstall:
	rm -f $(INCLUDEDIR)/adept.h
	rm -f $(LIBDIR)/libadept.a

clean:
	cd adept && $(MAKE) clean
	cd test && $(MAKE) clean
	cd advection && $(MAKE) clean
	cd multiscatter && $(MAKE) clean

clean-autosaves:
	rm -f *~
	rm -f */*~
	rm -f */*/*~

.PHONY: default all libadept.a multiscatter advection-all multiscatter-all install clean clean-autosaves