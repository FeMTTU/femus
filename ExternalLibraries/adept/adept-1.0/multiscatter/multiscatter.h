/* multiscatter.h -- Calculation of lidar and radar multiple scattering 

   Copyright (C) 2004-2011 Robin Hogan <r.j.hogan@reading.ac.uk>

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


   This is the header file for an algorithm for efficient calculation
   of the lidar or radar backscatter profile in the presence of
   multiple scattering.  Several different algorithms are available
   depending on the multiple scattering regime.
*/

#ifndef _MULTISCATTER_H
#define _MULTISCATTER_H 1

//#ifdef __cplusplus
//extern "C" {
//#endif

#include <stdio.h>

#define MS_VERSION "1.2"

/* *** PREPROCESSOR DEFINITIONS *** */

/* The code may be compiled either in single or double
   precision. Double precision is the default but may be overridden by
   compiling with "-DSINGLE_PRECISION". */
#ifdef SINGLE_PRECISION
typedef float ms_real;
#else
typedef double ms_real;
#endif

/* External functions return an integer with one of these values */
#define MS_SUCCESS 0
#define MS_FAILURE 1
#define MS_MEMORY_ALLOCATION_ERROR 2
#define MS_INPUT_ERROR 3
#define MS_INVALID_CONTEXT 4
#define MS_COMMAND_LINE_ERROR 5
#define MS_INPUT_FILE_ERROR 6

/* Bitwise options used in the "options" member of the ms_config
   structure */

/* General settings */
#define MS_QUIET (1L<<0)
#define MS_ANNULAR_DETECTORS (1L<<24)
#define MS_CHECK_FOR_NAN (1L<<25)
#define MS_NUMERICAL_JACOBIAN (1L<<8)

/* Settings for the MS_SMALL_ANGLE_PVC_ORIGINAL algorithm */
#define MS_SIMPLE_OPTICAL_DEPTH (1L<<2)
#define MS_CRUDE_OPTICAL_DEPTH (1L<<3)
#define MS_CRUDE_INTEGRATION (1L<<4)
#define MS_NO_MULTISCAT_WITHIN_GATE (1L<<5)
#define MS_DOUBLE_SCATTERING_ONLY (1L<<6)
#define MS_CRUDE_DOUBLE_SCATTERING (1L<<7)
#define MS_WIDE_ANGLE_CUTOFF (1L<<10)
#define MS_APPROXIMATE_EXPONENTIAL (1L<<11)
#define MS_OUTPUT_BSCATS (1L<<17)
#define MS_NO_MOLECULAR_EXTINCTION (1L<<9)

/* Settings for the MS_WIDE_ANGLE_TDTS algorithm */
#define MS_NO_FORWARD_LOBE (1L<<12)
#define MS_SSA_SCALES_FORWARD_LOBE (1L<<13)
#define MS_SIMPLE_2S_COEFFTS (1L<<14)
#define MS_PROPAGATION_TO_STDERR (1L<<15)
#define MS_IGNORE_SOURCE_WIDENING (1L<<26)

/* Settings for automatic differentiation */
#define MS_AUTOMATIC_DIFFERENTIATION_ADEPT (1L<<27)
#define MS_AUTOMATIC_DIFFERENTIATION_CPPAD (1L<<28)
#define MS_AUTOMATIC_DIFFERENTIATION_ADOLC (1L<<29)
#define MS_AUTOMATIC_DIFFERENTIATION_SACADO (1L<<30)
#define MS_AUTOMATIC_DIFFERENTIATION_SACADO_FAD (1L<<16)
#define MS_JACOBIAN (1L<<31)

/* Settings for the MS_SMALL_ANGLE_PVC_EXPLICIT algorithm */
#define MS_OUTPUT_DISTRIBUTION (1L<<18)

/* Specify the fraction by which the range spacing is allowed to vary
   before the data are regridded for the wide-angle algorithm. */
#define MS_RANGE_SPACING_TOLERANCE 0.05

/* The wavelength dividing the default radar and default lidar
   settings, in metres */
#define MS_RADAR_LIDAR_TRANSITION_WAVELENGTH 100.0e-6

/* *** STRUCTURES AND ENUMERATIONS *** */

/* The multiple scattering is split into two parts: (1)
   single-scattering + small-angle scattering, and (2) wide-angle
   scattering. Two enumerations are provided to select the possible
   algorithms that can be used for each part. */

/* First the single-scattering + small-angle scattering; note that PVC
   = "photon variance-covariance method", a framework for modelling
   the photon distribution. */
typedef enum {
  MS_SINGLE_AND_SMALL_ANGLE_NONE = 0, /* No single or small-angle scattering */
  MS_SINGLE_SCATTERING,        /* Single scattering only, e.g. for radar */
  MS_SMALL_ANGLE_PVC_ORIGINAL, /* Hogan (Applied Optics, 2006): speed O(N^2) */
  MS_SMALL_ANGLE_PVC_FAST,     /* Hogan (J Atmos Sci, 2008): O(N), DEFAULT */
  MS_SMALL_ANGLE_PVC_EXPLICIT, /* Eloranta-like but PVC framework: O(N^m) */
  MS_SMALL_ANGLE_PVC_FAST_LAG, /* As Hogan (2008) but calculating time
				  lag as well */
  MS_NUM_SMALL_ANGLE_ALGORITHMS
} ms_small_angle_algorithm;

/* Next the wide-angle scattering, where TDTS = "time-dependent
   two-stream method" */
typedef enum {
  MS_WIDE_ANGLE_NONE = 0, /* Now wide-angle scattering */
  MS_WIDE_ANGLE_TDTS_FORWARD_LOBE,   /* Hogan and Battaglia (J Atmos
					Sci, 2008), DEFAULT */
  MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE,/* Radar version of TDTS */
  MS_NUM_WIDE_ANGLE_ALGORITHMS
} ms_wide_angle_algorithm;

/* Configuration settings passed to functions */
typedef struct {
  ms_small_angle_algorithm small_angle_algorithm;
  ms_wide_angle_algorithm  wide_angle_algorithm;
  int options;
  int max_scattering_order; /* For explicit scattering orders */
  ms_real max_theta;        /* A small-angle option */
  int first_wide_angle_gate;/* 0=all,-1=auto,1+=specify */
  /* Michael Mishchenko argues that coherent backscatter should
     enhance the multiple scattering by a factor of 2 for radar, but
     not for lidar: up to the user to specify this */
  ms_real coherent_backscatter_enhancement;
  /* Outputs */
  ms_real *small_angle_lag; /* m */
  ms_real total_src;        /* Fraction of available energy injected
			       into wide-angle calculation */
  ms_real total_reflected;
  ms_real ss_multiplier;
  int force_jacobian; // -1: force reverse; 1: force forward; 0 default
} ms_config;

#define MS_DEFAULT_CONFIG {MS_SMALL_ANGLE_PVC_FAST, \
      MS_WIDE_ANGLE_TDTS_FORWARD_LOBE,		    \
      0, 4, 0.1, 0, 1.0, NULL, -1.0, -1.0, 1.0, 0}

/* Surface properties; note that at present these are not used, but
   the plan is to incorporate single and multiple scatterng from the
   surface */
typedef struct {
  ms_real sigma0;
  ms_real diffuse_albedo;
  ms_real direct_to_diffuse_albedo;
  ms_real diffuse_to_direct_backscatter;
  ms_real range; /* This might implicitly be at 0 or at the far end
		    of the ray... */
} ms_surface;

#define MS_DEFAULT_SURFACE {0.0, 0.0, 0.0, 0.0, 0.0}

/* Different receiver types */
typedef enum {
  MS_TOP_HAT = 0,
  MS_GAUSSIAN = 1
} ms_receiver_type;

/* Properties of the instrument */
typedef struct {
  ms_receiver_type receiver_type;
  ms_real altitude;   /* In m */
  ms_real wavelength; /* In m */
  ms_real rho_transmitter; /* transmitter_half-angle divergence;*/
  ms_real* rho_receiver;    /* reciever half-angle fields-of-view;*/
  int nfov; /* Number of fields-of-view (for multiple FOV lidar) */
} ms_instrument;

#define MS_DEFAULT_INSTRUMENT {MS_TOP_HAT, 0.0, 1.0, 1.0, NULL, 0}

/* For Fortran interface: store properties in one structure to be
   indexed by an integer */
typedef struct {
  ms_config config;
  ms_instrument instrument;
  ms_surface surface;
} ms_context;

#define MS_DEFAULT_CONTEXT {MS_DEFAULT_CONFIG, MS_DEFAULT_INSTRUMENT, MS_DEFAULT_SURFACE}

/* *** CONFIGURATION FUNCTIONS ***/

/* Functions to change the algorithm options */
void ms_set_options(ms_config* config, int options);
void ms_add_options(ms_config* config, int options);

#define MS_AUTO_FIRST_WIDE_ANGLE_GATE -1
void ms_set_first_wide_angle_gate(ms_config* config, int first_gate);

/* A version of Eloranta's (1998) small-angle multiple scattering
   algorithm but using the PVC methodology described by Hogan (2008),
   taken to "norder" orders of scattering */
int ms_set_max_scattering_order(ms_config* config, int norder);

/* Functions for setting and retrieving the photon spatial
   distribution from the MS_SMALL_ANGLE_EXPLICIT algorithm - note that
   these functions are not thread-safe */
int ms_set_distribution_resolution(int num_x, ms_real dx);
int ms_print_distribution(ms_config config, FILE* file);


/* MAIN INTERFACE FUNCTIONS */

/* Perform the multiple scattering calculation, using which ever
   combination of algorithms is specified in "config". If
   "bscat_air_out" is NULL then bscat_out will contain the sum of the
   particle and air backscatter. Otherwise the particle and air
   backscatters will be calculated separately and reported in
   bscat_out and bscat_air_out; this is useful for simulating
   measurements by Raman and HSRL lidars. */
int multiscatter(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,    /* Cloud/aerosol equivalent radius, microns */
    const ms_real *ext,       /* Total extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft., m-1 (NULL for vacuum) */
    const ms_real *ssa_air,   /* Air single-scatter albedo */
    const ms_real *droplet_fraction,/* Fraction of extinction from droplets */
    const ms_real *pristine_ice_fraction,/* Fraction of ext from pristine ice */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    ms_real *bscat_air_out    /* Measured backscatter of air, m-1 sr-1 */
    );

/* Multiple scattering calculation including returning the logarithmic
   Jacobian: USE WITH CARE; NO LONGER MAINTAINED; USE
   ms_jacobian_linear INSTEAD */
int ms_jacobian(
    /* Input data */
    int n,                    /* Number of range gates in profile */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    ms_real *range,           /* Height of each range gate, metres */
    ms_real *radius,          /* Cloud/aerosol equivalent radius, microns */
    ms_real *ext,             /* Cloud/aerosol extinction coefficient, m-1 */
    ms_real *ssa,             /* Total single-scatter albedo */
    ms_real *g,               /* Total asymmetry factor */
    ms_real *ext_bscat_ratio, /* Cloud/aerosol ext./backscatter ratio, sr */
    ms_real *ext_air,         /* Air ext. coefft, m-1 (NULL for vacuum) */
    ms_real *ssa_air,
    ms_real *bscat_peak_factor, /* Relative amplitude of P(180), 0-1 */
    ms_real *bscat_peak_width,  /* Width of peak near 180, radians */
    int n_x, /* Number of state variables to compute the Jacobian of */
    int n_y, /* Number of measurements to compute the Jacobian of */
    int *i_x, /* Index to the required state variables */
    int *i_y, /* Index to the required measurements */
    /* Output data */
    ms_real *bscat_out,        /* Measured backscatter, m-1 sr-1 */
    ms_real **d_ln_bscat_d_ln_ext, /* Jacobian with respect to ln(extinction) */
    ms_real **d_ln_bscat_d_ln_radius, /* Jacobian w.r.t. ln(radius) */
    ms_real *d_ln_bscat_d_ln_ext_bscat_ratio /* ...w.r.t. ln(ext_bscat_ratio) */
    );

/* Multiple scattering calculation including returning the linear
   Jacobian; in practice this uses a numerical technique to calcualate
   the Jacobian by perturbing each input variable in turn. It is
   therefore rather slow. The memory need to be allocated for the
   output Jacobians before calling this function (or use NULL for
   Jacobians that are not required). */
int ms_jacobian_linear(
    /* Input data */
    int n,                    /* Number of range gates in profile */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,    /* Cloud/aerosol equivalent radius, microns */
    const ms_real *ext,       /* Cloud/aerosol extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    const ms_real *ssa_air,
    const ms_real *droplet_fraction,/* Fraction of ext due to droplets, 0-1 */
    const ms_real *pristine_ice_fraction,/* Frac. of ext due to pristine ice */
    int n_x,  /* Number of state variables to compute the Jacobian of */
    int n_y,  /* Number of measurements to compute the Jacobian of */
    int *i_x, /* Index to the required state variables */
    int *i_y, /* Index to the required measurements */
    /* Output data */
    ms_real *bscat_out,            /* Measured backscatter, m-1 sr-1 */
    ms_real *bscat_air_out,        /* Measured backscatter, m-1 sr-1 */
    /* Output Jacobian matrices */
    ms_real **d_bscat_d_ext,       /* Jacobian with respect to ext */
    ms_real **d_bscat_d_ssa,       /* Jacobian w.r.t. radius */
    ms_real **d_bscat_d_g,         /* Jacobian with respect to ext */
    ms_real **d_bscat_d_radius,    /* Jacobian w.r.t. radius */
    ms_real *d_bscat_d_ext_bscat_ratio,/* ...w.r.t. ext_bscat_ratio */
    ms_real **d_bscat_air_d_ext,   /* Jacobian with respect to ext */
    ms_real **d_bscat_air_d_ssa,   /* Jacobian w.r.t. radius */
    ms_real **d_bscat_air_d_g,     /* Jacobian with respect to ext */
    ms_real **d_bscat_air_d_radius /* Jacobian w.r.t. radius */
 );


/* Adjoint of the multiple scattering algorithm */
int multiscatter_AD(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,    /* Cloud/aerosol equivalent radius, microns */
    const ms_real *ext,       /* Total extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft., m-1 (NULL for vacuum) */
    const ms_real *ssa_air,   /* Air single-scatter albedo */
    const ms_real *droplet_fraction, /* Fraction of extinction from droplets */
    const ms_real *pristine_ice_fraction,/* Fraction of ext from pristine ice */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    ms_real *bscat_air_out,   /* Measured backscatter of air, m-1 sr-1 */
    /* Adjoint inputs */
    const ms_real *bscat_AD,
    const ms_real *bscat_air_AD,
    /* Adjoint outputs */
    ms_real *radius_AD,
    ms_real *ext_AD,
    ms_real *ssa_AD,
    ms_real *g_AD,
    ms_real *ext_bscat_ratio_AD,
    ms_real *ext_air_AD,
    ms_real *droplet_fraction_AD,
    ms_real *pristine_ice_fraction_AD);


/* *** DIAGNOSTIC FUNCTIONS *** */

/* Return 1 if any of the n elements of x contains "not a number", 0
   otherwise.  This is useful to check the inputs or outputs of the
   algorithm. */
int ms_isnan(int n, const ms_real* x);

/* Write an "input" file from the variables provided; this is useful
   when the multiscatter algorithm is embedded within a larger
   retrieval algorithm, a particular profile has produced NaN outputs
   and it is necessary to know what the inputs were. */
int ms_dump_input(
    const char* filename, /* Name of file to save data to */
    /* Input data */
    int n,                    /* Number of input gates */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,    /* Cloud/aerosol equivalent radius, microns */
    const ms_real *ext,       /* Total extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft., m-1 (NULL for vacuum) */
    const ms_real *ssa_air,   /* Air single-scatter albedo */
    const ms_real *droplet_fraction,/* Fraction of extinction from droplets */
    const ms_real *pristine_ice_fraction,/* Fraction of ext from pristine ice */
    /* Adjoint inputs */
    const ms_real *bscat_AD,
    const ms_real *bscat_air_AD);

/* Print the characteristics of the algorithms that will be run to the
   specified file (e.g. stderr) */
void ms_print_algorithms(ms_config config,
			 ms_instrument instrument,
			 const ms_real* range, 
			 int use_isotropic_pp,
			 int separate_bscat_air,
			 int output_adjoint,
			 int output_jacobian,
			 int calc_jacobian,
			 FILE* file);


/* Check that range gate spacing is monotonic and that the spacings
   are the same to within the specified tolerance (e.g. 0.05 for 5% or
   0.01 for 1%). Return 1 if regular, 0 if irregular */
int ms_range_spacing_is_regular(int n, const ms_real* range, 
				ms_real tolerance);

/* Prints statistics to the specified stream */
int ms_print_stats(ms_config config, FILE *file);

int
simple_algorithm(
    /* Input data */
    int n,                   
    int m,                   
    ms_config *config,       
    ms_instrument instrument,
    ms_surface surface,   
    const ms_real *dummy1,
    const ms_real *a,
    const ms_real *b,
    const ms_real *c,
    const ms_real *d,
    const ms_real *e,
    const ms_real *dummy2,
    const ms_real *dummy3,
    const ms_real *f,
    const ms_real *g,
    /* Output data */
    ms_real *out,      
    ms_real *dummy_out);


//#ifdef __cplusplus
//}
//#endif

#endif
