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

#ifndef _MS_AD_H
#define _MS_AD_H 1

#include <stdio.h>

#include "multiscatter.h"

/* Simple single scattering calculation with no delta-Eddington
   scaling, usually used in conjunction with the TDTS multiple
   scattering algorithm */
int ms_singlescatter(
    /* Input data */
    int n,                    /* Number of range gates in profile */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const areal *ext,       /* Cloud/aerosol extinction coefficient, m-1 */
    const areal *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    /* Output data */
    areal *bscat_out,       /* Measured backscatter, m-1 sr-1 */ 
    areal *bscat_air_out);

/* The Hogan (2008) fast small-angle multiple scattering algorithm is
   called using this function. It is O(n) efficient. MS_FAILURE is
   returned if there is an error, MS_SUCCESS otherwise. */
int ms_fast_small_angle(
    /* Input data */
    int n,                    /* Number of range gates in profile */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const areal *radius,    /* Cloud/aerosol equivalent radius, microns */
    const areal *ext,       /* Cloud/aerosol extinction coefficient, m-1 */
    const areal *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    const areal *droplet_fraction,/* Fraction of ext due to droplets, 0-1 */
    const areal *pristine_ice_fraction, /* ...due to pristine ice with
					     Yang-like phase function */
    /* Output data */
    areal *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    areal *bscat_air_out);

/* A version of Eloranta's (1998) small-angle multiple scattering
   algorithm but using the PVC methodology described by Hogan (2008),
   taken to "norder" orders of scattering */
int ms_explicit(
    int n,                   /* Number of range gates in profile */
    ms_config *config,       /* Configuration information */
    ms_instrument instrument,/* Structure containing instrument variables */
    ms_surface surface,      /* Surface scattering variables */
    const ms_real *range,    /* Height of each range gate, metres */
    const areal *radius,   /* Cloud/aerosol equivalent radius, microns */
    const areal *ext,      /* Cloud/aerosol extinction coefficient, m-1 */
    const areal *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,  /* Air ext. coefft, m-1 (NULL for vacuum) */
    areal *bscat_out,      /* Measured backscatter, m-1 sr-1 */
    areal *bscat_air_out);


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
    const areal *ext,       /* Total extinction coefficient, m-1 */
    const areal *ssa,       /* Total single-scatter albedo */
    const areal *g,         /* Total asymmetry factor */
    const areal *src_power_in,/* Source function power: inwards */
    const areal *src_power_out,/* Source function power: outwards */
    const areal *src_width2,/* 1/e source function width^2, radians^2 */
    /* Output data */
    areal *bscat_out);      /* Measured backscatter, m-1 sr-1 */

/* Just do the wide-angle part of the calculation */
int ms_wide_angle(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const areal *radius,    /* Height of each range gate, metres */
    const areal *ext,       /* Total extinction coefficient, m-1 */
    const areal *ssa,       /* Total single-scatter albedo */
    const areal *g,         /* Total asymmetry factor */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    const ms_real *ssa_air,
    /* Output data */
    areal *bscat_out);      /* Measured backscatter, m-1 sr-1 */

int ms_wide_angle_regrid(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const areal *radius,     /* Height of each range gate, metres */
    const areal *ext,       /* Total extinction coefficient, m-1 */
    const areal *ssa,       /* Total single-scatter albedo */
    const areal *g,         /* Total asymmetry factor */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    const ms_real *ssa_air,
    /* Output data */
    areal *bscat_out);       /* Measured backscatter, m-1 sr-1 */

/* Calculate and return the scaling factor to account for anisotropic
   near-backscatter phase functions */
areal ms_anisotropic_factor(
   areal width2,          /* spatial variance (m2) */
   areal zeta2,           /* angular variance (rad2) */
   areal cov,             /* covariance (m rad) */
   areal Theta2,          /* forward lobe variance (rad2) */
   ms_real range,           /* range from instrument (m) */
   ms_real width2_max,      /* receiver field-of-view variance (m2) */
   areal ext,             /* Cloud extinction coefficient (m-1) */
   ms_real ext_air,         /* Air extinction coefficient (m-1) */
   areal droplet_fraction,/* Fraction of ext due to droplets */
   areal pristine_ice_fraction);/* Fraction of ext due to pristine ice */


/* Calculate the spatial variance of the quasi-direct outgoing beam of
   photons using an O(N) efficient algorithm */
int ms_variance(int n, ms_real wavelength, ms_real rho_laser,
		ms_real lidar_altitude, const ms_real *range,
		const areal *radius, const areal *ext,
		areal *variance_out);

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
    const areal *radius,    /* Cloud/aerosol equivalent radius, microns */
    const areal *ext,       /* Total extinction coefficient, m-1 */
    const areal *ssa,       /* Total single-scatter albedo */
    const areal *g,         /* Total asymmetry factor */
    const areal *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft., m-1 (NULL for vacuum) */
    const ms_real *ssa_air,   /* Air single-scatter albedo */
    const areal *droplet_fraction,/* Fraction of extinction from droplets */
    const areal *pristine_ice_fraction,/* Fraction of ext from pristine ice */
    /* Output data */
    areal *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    areal *bscat_air_out    /* Measured backscatter of air, m-1 sr-1 */
    );


int
simple_algorithm(
    /* Input data */
    int n,                   
    int m,                   
    ms_config *config,       
    ms_instrument instrument,
    ms_surface surface,   
    const ms_real *dummy1,
    const areal *a,
    const areal *b,
    const areal *c,
    const areal *d,
    const areal *e,
    const ms_real *dummy2,
    const ms_real *dummy3,
    const areal *f,
    const areal *g,
    /* Output data */
    areal *out,      
    areal *dummy_out);



#endif
