// advection_schemes_AD.h - Header for the hand-coded adjoints

// Hand-coded adjoint of Lax-Wendroff advection scheme
void lax_wendroff_AD(int nt, double c, const double q_init[NX], double q[NX],
		     const double q_AD_const[NX], double q_init_AD[NX]);

// Hand-coded adjoint of Toon advection scheme
void toon_AD(int nt, double c, const double q_init[NX], double q_out[NX],
	     const double q_AD_const[NX], double q_init_AD[NX]);
