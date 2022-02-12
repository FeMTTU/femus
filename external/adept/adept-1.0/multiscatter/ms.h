/* ms.h -- Calculation of lidar and radar multiple scattering 

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


   This is the header file for the internal functions of the
   multiscatter library; external programs should only need to include
   multiscatter.h.
*/

#ifndef _MS_H
#define _MS_H 1

// Need to allow C++ overloading
//#ifdef __cplusplus
//extern "C" {
//#endif

#include <stdio.h>

#include "multiscatter.h"

#define MS_PI 3.14159265358979323846

/* SETTINGS FOR TIME-DEPENDENT TWO-STREAM ALGORITHM */

/* This is from diffusion theory and should not be changed */
#define MS_D_FACTOR 1.33333333
#define MS_ONE_OVER_D_FACTOR 0.75

/* The speed of light */
#define MS_C (2.99792458e8)

/* The two-stream cosine-zenith angle */
#define MS_MU1 0.5

/* For speed it is useful not to calculate multiple scattering on the
   profile, but only the part that is optically thick enough that
   multiple scattering will be significant. The threshold is in terms
   of the ratio of a mean-free-path to the radar/lidar footprint - if
   it is less than this value then multiple scattering will be
   calculated for all gates beyond this point. A typical number would
   be 20, but if it is set very large then multiple scattering will be
   calculated for the entire profile. */
#define MS_MFP_FOOTPRINT_THRESHOLD 50.0

/* The input data may be on an irregular grid, while the multiple
   scattering calculation is performed on a regular grid. In deciding
   the resolution to use for this regular grid, there is a trade-off
   between the accuracy and the speed. However, there is no need to
   run the algorithm at a higher resolution than the important initial
   region of the cloud, and after many optical depths into the cloud
   all sharp gradients are smoothed out anyway. To decide what
   resolution to use, the cloud is searched up to a particular
   scattering optical depth, and the finest grid spacing in this
   region is used for the regular grid. This variable determines the
   optical depth threshold to use. A typical value would be 5. */
#define MS_OD_DRANGE_THRESHOLD 5.0

/* SETTINGS FOR PHOTON VARIANCE-COVARIANCE METHOD */

/* Factors parameterizing the shape of the droplet phase function near
   backscatter */
#define MS_MIE_U1 0.3
#define MS_MIE_U2 0.5
#define MS_MIE_V1_SQD 16.0
#define MS_MIE_V2_SQD 0.16
/* Factors parameterizing the shape of the droplet phase function near
   backscatter */
#define MS_YANG_W 0.89
#define MS_YANG_GAMMA0 0.038

/* After calling ms_small_angle(), the following variables are
   available to be accessed - see multiscatter_qsa.c for details. */
extern ms_real *ms_E_once;           /* dimensionless */
extern ms_real *ms_Emu2_once;        /* radians^2 */
extern ms_real *ms_E_multi;          /* dimensionless */
extern ms_real *ms_Emu2_multi;       /* radians^2 */
extern ms_real *ms_Pmu2_once;        /* radians^2 */
extern ms_real *ms_Pmu2_multi;       /* radians^2 */
extern ms_real *ms_EPcov_once;       /* radians^2 */
extern ms_real *ms_EPcov_multi;      /* radians^2 */
extern ms_real *ms_bscat_single;     /* m-1 sr-1 */
extern ms_real *ms_bscat_double;     /* m-1 sr-1 */
extern ms_real *ms_bscat_multi;      /* m-1 sr-1 */
extern ms_real *ms_bscat_air_single; /* m-1 sr-1 */
extern ms_real *ms_bscat_air_double; /* m-1 sr-1 */
extern ms_real *ms_bscat_air_multi;  /* m-1 sr-1 */


/* Functions to allocate and zero or to free intermediate
   arrays. External programs need not call these functions.*/
int ms_init_intermediate_arrays(int ngates);
void ms_free_intermediate_arrays();


/* Simple single scattering calculation with no delta-Eddington
   scaling, usually used in conjunction with the TDTS multiple
   scattering algorithm */
int ms_singlescatter(
    /* Input data */
    int n,                    /* Number of range gates in profile */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *ext,       /* Cloud/aerosol extinction coefficient, m-1 */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */ 
    ms_real *bscat_air_out);

/* The Hogan (2006) small-angle multiple scattering algorithm is
   called using this function. It is O(n^2) efficient. MS_FAILURE is
   returned if there is a error, MS_SUCCESS otherwise. */
int ms_small_angle(
    /* Input data */
    int n,                    /* Number of range gates in profile */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,    /* Cloud/aerosol equivalent radius, microns */
    const ms_real *ext,       /* Cloud/aerosol extinction coefficient, m-1 */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft., m-1 (NULL for vacuum) */
    const ms_real *droplet_fraction, /* NOT YET IMPLEMENTED */
    const ms_real *pristine_ice_fraction,  /* NOT YET IMPLEMENTED */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    ms_real *bscat_air_out);    /* Measured backscatter of air, m-1 sr-1 */

/* The Hogan (2008) fast small-angle multiple scattering algorithm is
   called using this function. It is O(n) efficient. MS_FAILURE is
   returned if there is an error, MS_SUCCESS otherwise. */
int ms_fast_small_angle(
    /* Input data */
    int n,                    /* Number of range gates in profile */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,    /* Cloud/aerosol equivalent radius, microns */
    const ms_real *ext,       /* Cloud/aerosol extinction coefficient, m-1 */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    const ms_real *droplet_fraction,/* Fraction of ext due to droplets, 0-1 */
    const ms_real *pristine_ice_fraction, /* ...due to pristine ice with
					     Yang-like phase function */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    ms_real *bscat_air_out);

int ms_fast_small_angle_AD(
    /* Input data */
    int n,                    /* Number of range gates in profile */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,    /* Cloud/aerosol equivalent radius, microns */
    const ms_real *ext,       /* Cloud/aerosol extinction coefficient, m-1 */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    const ms_real *droplet_fraction,/* Fraction of ext due to droplets, 0-1 */
    const ms_real *pristine_ice_fraction, /* ...due to pristine ice with
					     Yang-like phase function */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    ms_real *bscat_air_out,   /* Measured backscatter of air, m-1 sr-1 */
    /* Adjoint inputs */
    const ms_real *bscat_AD,
    const ms_real *bscat_air_AD,
    /* Adjoint outputs */
    ms_real *radius_AD,
    ms_real *ext_AD,
    ms_real *ext_bscat_ratio_AD,
    ms_real *ext_air_AD,
    ms_real *droplet_fraction_AD,
    ms_real *pristine_ice_fraction_AD);

int
ms_fast_small_angle_extras(
    /* Input data */
    int n,                    /* Number of range gates in profile */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    ms_real *range,           /* Height of each range gate, metres */
    ms_real *radius,          /* Cloud/aerosol equivalent radius, microns */
    ms_real *ext,             /* Cloud/aerosol extinction coefficient, m-1 */
    ms_real *ext_bscat_ratio, /* Cloud/aerosol ext./backscatter ratio, sr */
    ms_real *ext_air,         /* Air ext. coefft, m-1 (NULL for vacuum) */
    ms_real *droplet_fraction,/* Fraction of ext due to droplets, 0-1 */
    ms_real *pristine_ice_fraction, /* ...due to pristine ice with Yang-like
				 phase function */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    ms_real *bscat_air_out,   /* Measured backscatter of air, m-1 sr-1 */
 /* Transmission Enhancement factor due to multiple scattering but
    neglecting anisotropic phase functions: this is required for the
    approximate adjoint (can be NULL): */
    ms_real *trans_enhancement_out,
    ms_real *anisotropic_factor_out); /* Anisotropic backscatter factor */

/* As ms_fast_small_angle but with a pulse lag
   calculation */
int ms_fast_small_angle_lag(
    /* Input data */
    int n,                    /* Number of range gates in profile */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,    /* Cloud/aerosol equivalent radius, microns */
    const ms_real *ext,       /* Cloud/aerosol extinction coefficient, m-1 */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    const ms_real *droplet_fraction,/* Fraction of ext due to droplets, 0-1 */
    const ms_real *pristine_ice_fraction, /* ...due to pristine ice with
					     Yang-like phase function */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    ms_real *bscat_air_out,   /* Measured backscatter of air, m-1 sr-1 */
    ms_real *lag_out);        /* Pulse lag, m */

/* A version of Eloranta's (1998) small-angle multiple scattering
   algorithm but using the PVC methodology described by Hogan (2008),
   taken to "norder" orders of scattering */
int ms_explicit(
    int n,                   /* Number of range gates in profile */
    ms_config *config,       /* Configuration information */
    ms_instrument instrument,/* Structure containing instrument variables */
    ms_surface surface,      /* Surface scattering variables */
    const ms_real *range,    /* Height of each range gate, metres */
    const ms_real *radius,   /* Cloud/aerosol equivalent radius, microns */
    const ms_real *ext,      /* Cloud/aerosol extinction coefficient, m-1 */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,  /* Air ext. coefft, m-1 (NULL for vacuum) */
    ms_real *bscat_out,      /* Measured backscatter, m-1 sr-1 */
    ms_real *bscat_air_out);


/* Perform the Hogan and Battaglia (2008) time-dependent two-stream
   calculation for wide-angle multiple scattering, and add the result
   to bscat_out */
int ms_tdts(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *ext,       /* Total extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *src_power_in,/* Source function power: inwards */
    const ms_real *src_power_out,/* Source function power: outwards */
    const ms_real *src_width2,/* 1/e source function width^2, radians^2 */
    /* Output data */
    ms_real *bscat_out);      /* Measured backscatter, m-1 sr-1 */

/* Just do the wide-angle part of the calculation */
int ms_wide_angle(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,    /* Height of each range gate, metres */
    const ms_real *ext,       /* Total extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    const ms_real *ssa_air,
    /* Output data */
    ms_real *bscat_out);      /* Measured backscatter, m-1 sr-1 */

int ms_wide_angle_regrid(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,     /* Height of each range gate, metres */
    const ms_real *ext,       /* Total extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    const ms_real *ssa_air,
    /* Output data */
    ms_real *bscat_out);       /* Measured backscatter, m-1 sr-1 */

/* Calculate and return the scaling factor to account for anisotropic
   near-backscatter phase functions */
ms_real ms_anisotropic_factor(
   ms_real width2,          /* spatial variance (m2) */
   ms_real zeta2,           /* angular variance (rad2) */
   ms_real cov,             /* covariance (m rad) */
   ms_real Theta2,          /* forward lobe variance (rad2) */
   ms_real range,           /* range from instrument (m) */
   ms_real width2_max,      /* receiver field-of-view variance (m2) */
   ms_real ext,             /* Cloud extinction coefficient (m-1) */
   ms_real ext_air,         /* Air extinction coefficient (m-1) */
   ms_real droplet_fraction,/* Fraction of ext due to droplets */
   ms_real pristine_ice_fraction);/* Fraction of ext due to pristine ice */

ms_real ms_anisotropic_factor_AD(
   ms_real width2,          /* spatial variance (m2) */
   ms_real zeta2,           /* angular variance (rad2) */
   ms_real cov,             /* covariance (m rad) */
   ms_real Theta2,          /* forward lobe variance (rad2) */
   ms_real range,           /* range from instrument (m) */
   ms_real width2_max,      /* receiver FOV variance (m2) */
   ms_real bscat,           /* Cloud bscat coefficient (m-1) */
   ms_real bscat_air,       /* Air bscat coefficient (m-1) */
   ms_real droplet_fraction,/* Fraction of ext from droplets */
   ms_real pristine_ice_fraction, /* Fraction of ext from pristine ice */
   /* Adjoint input */
   ms_real factor_AD,       /* Adjoint of anisotropic factor */
   /* Adjoint outputs */
   ms_real* width2_AD,
   ms_real* zeta2_AD,
   ms_real* cov_AD,
   ms_real* Theta2_AD,
   ms_real* bscat_AD,
   ms_real* bscat_air_AD,
   ms_real* droplet_fraction_AD,
   ms_real* pristine_ice_fraction_AD);



/* Calculate the spatial variance of the quasi-direct outgoing beam of
   photons using an O(N) efficient algorithm */
int ms_variance(int n, ms_real wavelength, ms_real rho_laser,
		ms_real lidar_altitude, const ms_real *range,
		const ms_real *radius, const ms_real *ext,
		ms_real *variance_out);

/* Functions for setting and retrieving the photon spatial
   distribution - note that these functions are not thread-safe */
int ms_init_distribution_arrays(int num_gates);
int ms_increment_distribution(int order, int gate, 
			       ms_real energy, ms_real variance);


/* ADJOINT FUNCTIONS */

int ms_wide_angle_AD(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,     /* Height of each range gate, metres */
    const ms_real *ext,       /* Total extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    const ms_real *ssa_air,
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    /* Adjoint inputs */
    const ms_real *bscat_AD,
    /* Adjoint outputs */
    ms_real *radius_AD,
    ms_real *ext_AD,
    ms_real *ssa_AD,
    ms_real *g_AD);

int ms_tdts_AD(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *ext,       /* Total extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *src_power_in,/* Source function power: inwards */
    const ms_real *src_power_out,/* Source function power: outwards */
    const ms_real *src_width2,/* 1/e source function width^2, radians^2 */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    /* Adjoint inputs */
    const ms_real *bscat_AD,
    /* Adjoint outputs */
    ms_real *ext_AD,
    ms_real *ssa_AD,
    ms_real *g_AD,
    ms_real *src_power_in_AD,
    ms_real *src_power_out_AD,
    ms_real *src_width2_AD);

int ms_singlescatter_AD(
    /* Input data */
    int n,                    /* Number of range gates in profile */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *ext,       /* Cloud/aerosol extinction coefficient, m-1 */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */ 
    ms_real *bscat_air_out,
    /* Adjoint terms (in) */
    const ms_real *bscat_AD,
    const ms_real *bscat_air_AD,
    /* Adjoint terms (out) */
    ms_real *ext_AD,
    ms_real *ext_bscat_ratio_AD);

int ms_small_angle_ADonly(int n, ms_instrument instrument,
		    ms_surface surface, ms_real *range,
		    ms_real *ext, ms_real *ext_bscat_ratio,
		    ms_real *ext_air, ms_real *ssa_air,
		    ms_real *bscat, ms_real *bscat_air_out,
		    ms_real *trans_enhancement,
		    ms_real *anisotropic_factor,
		    /* Adjoint terms (in) */
		    ms_real *bscat_AD, ms_real *bscat_air_AD,
		    /* Adjoint terms (out) */
		    ms_real *ext_AD, ms_real *ext_bscat_ratio_AD);

int
ms_variance_AD(int n, ms_real wavelength, ms_real rho_laser,
	       ms_real inst_altitude, const ms_real *range, 
	       const ms_real *radius, const ms_real *ext, 
	       ms_real *variance_out,
	       ms_real *variance_AD,
	       ms_real *radius_AD, ms_real *ext_AD);

ms_real ms_get_drange(int i, int n, const ms_real *range);
ms_real ms_get_midpoint(int i, int n, const ms_real *range);



//#ifdef __cplusplus
//}
//#endif

#endif
