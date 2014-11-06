// test_gsl_interface.cpp - "main" function for Test 4

// This program minimizes the N-dimensional Rosenbrock banana
// function, with the number of dimensions optionally provided on the
// command line

#include <iostream>
#include <vector>

#include "state.h"

int
main(int argc, char** argv)
{
  std::cout << "Testing Adept-GSL interface using N-dimensional Rosenbrock function\n";
  std::cout << "Usage: " << argv[0] << " [number_of_dimensions]\n";

  // Read number of dimensions from the command line (default 2)
  int nx = 2;
  if (argc > 1) {
    nx = atoi(argv[1]);
  }
   
  if (nx < 2) {
    std::cout << "Error: must have 2 or more dimensions, but "
	      << nx << " requested\n";
    exit(1);
  }

  // Create minimization environment (see state.h) and then minimize
  // the function; note that initial values are set on construction.
  State state(nx);
  state.minimize();

  // Print out the result
  std::vector<double> x;
  state.x(x);
  std::cout << "Final state: x = [";
  for (unsigned int i = 0; i < nx; i++) {
    std::cout << " " << x[i];
  }
  std::cout << "]\n";
  
  return 0;
}
