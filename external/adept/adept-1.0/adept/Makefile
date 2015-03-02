# Most of the library is in the adept.h header file since it needs to
# be inlined for efficiency, but the rest of it is in adept.cpp, which
# is compiled here into a static library.

include ../Makefile_include
OBJECTS = adept.o
LIB = ../lib/libadept.a
#VECTORIZE_FLAGS = -ftree-vectorizer-verbose=5

all: $(LIB)

$(LIB): $(OBJECTS)
	ar r $(LIB) $(OBJECTS)

%.o: %.cpp ../include/adept*.h
	$(CPP) $(CPPFLAGS) $(VECTORIZE_FLAGS) $(INCLUDES) -c $<

clean:
	rm -f $(OBJECTS) $(LIB)

clean-autosaves:
	rm -f *~

