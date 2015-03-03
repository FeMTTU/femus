// advection.h - Header file for two advection algorithms from the Adept paper

#define NX 100

// Lax-Wendroff advection scheme
void lax_wendroff(int nt, double courant, const adouble q_init[NX], adouble q[NX]);

// Toon advection scheme
void toon(int nt, double courant, const adouble q_init[NX], adouble q[NX]);

