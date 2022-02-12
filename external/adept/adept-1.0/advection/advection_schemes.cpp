// advection.cpp - Two test advection algorithms from the Adept paper

// This code gets compiled multiple times with different automatic
// differentiation tools in order that their performance can be
// compared.

#include <iostream>
#include <iomanip>
#include <cmath>
#include "ad_library_switch.h"

//using namespace adept;

// Lax-Wendroff scheme applied to linear advection
void lax_wendroff(int nt, double c, const adouble q_init[NX], adouble q[NX]) {
  preallocate_statements((nt+1)*NX*3);
  preallocate_operations((nt+1)*NX*7);
  adouble flux[NX-1];                        // Fluxes between boxes
  for (int i=0; i<NX; i++) q[i] = q_init[i]; // Initialize q 
  for (int j=0; j<nt; j++) {                 // Main loop in time
    for (int i=0; i<NX-1; i++) flux[i] = 0.5*c*(q[i]+q[i+1]+c*(q[i]-q[i+1]));
    for (int i=1; i<NX-1; i++) q[i] += flux[i-1]-flux[i];
    q[0] = q[NX-2]; q[NX-1] = q[1];          // Treat boundary conditions
  }
}

// Toon advection scheme applied to linear advection
void toon(int nt, double c, const adouble q_init[NX], adouble q[NX]) {
  preallocate_statements((nt+1)*NX*3);
  preallocate_operations((nt+1)*NX*9);
  adouble flux[NX-1];                        // Fluxes between boxes
  for (int i=0; i<NX; i++) q[i] = q_init[i]; // Initialize q
  for (int j=0; j<nt; j++) {                 // Main loop in time
    for (int i=0; i<NX-1; i++) flux[i] = (exp(c*log(q[i]/q[i+1]))-1.0) 
                                         * q[i]*q[i+1] / (q[i]-q[i+1]);
    for (int i=1; i<NX-1; i++) q[i] += flux[i-1]-flux[i];
    q[0] = q[NX-2]; q[NX-1] = q[1];          // Treat boundary conditions
  }
}
