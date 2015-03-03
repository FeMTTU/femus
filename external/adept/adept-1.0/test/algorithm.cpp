// algorithm.cpp - A simple demonstration algorithm used in Tests 1 & 2 

#include <cmath>

#include "algorithm.h"
using adept::adouble;

// A simple demonstration algorithm used in the Adept paper. Note that
// this algorithm can be compiled with
// -DADEPT_NO_AUTOMATIC_DIFFERENTIATION to create a version that takes
// double arguments and returns a double result.
adouble algorithm(const adouble x[2]) {
  adouble y = 4.0;
  adouble s = 2.0*x[0] + 3.0*x[1]*x[1];
  y *= sin(s);
  return y;
}
 
