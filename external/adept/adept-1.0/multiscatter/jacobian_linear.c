/* ms_jacobian_linear.c -- Linear numerical Jacobian of multiscatter algorithm

   Copyright (C) 2005-2010 Robin Hogan <r.j.hogan@reading.ac.uk>

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

*/

#include <stdlib.h>
#include <math.h>

#include "multiscatter.h"

/* Calculate the Jacobian of the multiscatter algorithm by rerunning
   the algorithm for perturbations of the input values. This is slow
   but reliable. */
static
int
ms_numerical_jacobian_linear(
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
    const ms_real *pristine_ice_fraction,/* Fraction of ext due to pristine ice */
    int n_x,  /* Number of state variables to compute the Jacobian of */
    int n_y,  /* Number of measurements to compute the Jacobian of */
    int *i_x, /* Index to the required state variables */
    int *i_y, /* Index to the required measurements */
    /* Output data: each Jacobian may be NULL if not required */
    ms_real *bscat,               /* Measured backscatter, m-1 sr-1 */
    ms_real *bscat_air,
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
    )
{
  ms_real tweak = 10000; /* Around 0.01% */
  ms_real exp_tweak = exp(1.0/tweak);
  ms_real bscat_tweak[m*instrument.nfov];
  ms_real bscat_air_tweak_data[m*instrument.nfov];
  ms_real* bscat_air_tweak = NULL;
  ms_real x[n];
  int i;
  /* Ordinary run */
  multiscatter(n, m, config, instrument, surface, range, radius,
	       ext, ssa, g, ext_bscat_ratio, ext_air, ssa_air,
	       droplet_fraction, pristine_ice_fraction,
	       bscat, bscat_air);


  /* Set each element of each array to zero. */
  for (i = 0; i < m*instrument.nfov; i++) {
    bscat_tweak[i] = 0.0;
  }

  /* If bscat_air is non-NULL then we need to store the output of the
     molecular channel when the inputs are tweaked */
  if (bscat_air) {
    bscat_air_tweak = bscat_air_tweak_data;
    for (i = 0; i < m*instrument.nfov; i++) {
      bscat_air_tweak_data[i] = 0.0;
    }
  }

  /* Set the tweaked inputs to zero */
  for (i = 0; i < n; i++) {
    x[i] = 1.0;
  }


  /* JACOBIAN WITH RESPECT TO EXTINCTION-TO-BACKSCATTER RATIO */

  if (d_bscat_d_ext_bscat_ratio) {
    /* Copy ext_bscat_ratio into x */
    int i, k, l;
    for (i = 0; i < n; i++) {
      x[i] = ext_bscat_ratio[i];
    }

    for (l = 0; l < n_y; l++) {
      d_bscat_d_ext_bscat_ratio[l] = 0.0;
    }
    
    /* Loop through all the x values to find the derivative for */
    for (k = 0; k < n_x; k++) {
      x[i_x[k]] *= exp_tweak;

      multiscatter(n, m, config, instrument, surface, range, radius,
		   ext, ssa, g, x, ext_air, ssa_air,
		   droplet_fraction, pristine_ice_fraction,
		   bscat_tweak, bscat_air_tweak);
      /* Loop through all the required measurements, which may be in
	 the multiple fields-of-view */
      for (l = 0; l < n_y; l++) {
	d_bscat_d_ext_bscat_ratio[l] += (bscat_tweak[i_y[l]] - bscat[i_y[l]])
	  / (x[i_x[k]] - ext_bscat_ratio[i_x[k]]);
      }
      x[i_x[k]] = ext_bscat_ratio[i_x[k]];
    }
  }


  /* JACOBIAN WITH RESPECT TO EXTINCTION */

  if (d_bscat_d_ext) {
    /* Copy ext into x */
    int i, k, l;
    for (i = 0; i < n; i++) {
      x[i] = ext[i];
    }

    /* Loop through all the x values to find the derivative for */
    for (k = 0; k < n_x; k++) {
      if (x[i_x[k]] != 0.0) {
	x[i_x[k]] *= exp_tweak;
      }
      else {
	x[i_x[k]] += 1.0/tweak;
      }

      multiscatter(n, m, config, instrument, surface, range, radius,
		   x, ssa, g, ext_bscat_ratio, ext_air, ssa_air,
		   droplet_fraction, pristine_ice_fraction,
		   bscat_tweak, bscat_air_tweak);
      /* Loop through all the required measurements, which may be in
	 the multiple fields-of-view */
      for (l = 0; l < n_y; l++) {
	d_bscat_d_ext[k][l] = (bscat_tweak[i_y[l]] - bscat[i_y[l]])
	  / (x[i_x[k]] - ext[i_x[k]]);
	if (bscat_air) {
	  d_bscat_air_d_ext[k][l] 
	    = (bscat_air_tweak[i_y[l]] - bscat_air[i_y[l]])
	    / (x[i_x[k]] - ext[i_x[k]]);
	}
      }
      x[i_x[k]] = ext[i_x[k]];
    }
  }



  /* JACOBIAN WITH RESPECT TO SINGLE-SCATTERING ALBEDO */

  if (d_bscat_d_ssa) {
    /* Copy ssa into x */
    int i, k, l;
    for (i = 0; i < n; i++) {
      x[i] = ssa[i];
    }

    /* Loop through all the x values to find the derivative for */
    for (k = 0; k < n_x; k++) {
      if (x[i_x[k]] != 0.0) {
	x[i_x[k]] /= exp_tweak;
      }
      else {
	x[i_x[k]] += 1.0/tweak;
      }

      multiscatter(n, m, config, instrument, surface, range, radius,
		   ext, x, g, ext_bscat_ratio, ext_air, ssa_air,
		   droplet_fraction, pristine_ice_fraction,
		   bscat_tweak, bscat_air_tweak);
      /* Loop through all the required measurements, which may be in
	 the multiple fields-of-view */
      for (l = 0; l < n_y; l++) {
	d_bscat_d_ssa[k][l] = (bscat_tweak[i_y[l]] - bscat[i_y[l]])
	  / (x[i_x[k]] - ssa[i_x[k]]);
	if (bscat_air) {
	  d_bscat_air_d_ssa[k][l] 
	    = (bscat_air_tweak[i_y[l]] - bscat_air[i_y[l]])
	    / (x[i_x[k]] - ssa[i_x[k]]);
	}
      }
      x[i_x[k]] = ssa[i_x[k]];
    }
  }


  /* JACOBIAN WITH RESPECT TO ASYMMETRY FACTOR */

  if (d_bscat_d_g) {
    /* Copy g into x */
    int i, k, l;
    for (i = 0; i < n; i++) {
      x[i] = g[i];
    }

    /* Loop through all the x values to find the derivative for */
    for (k = 0; k < n_x; k++) {
      if (x[i_x[k]] != 0.0) {
	x[i_x[k]] *= exp_tweak;
      }
      else {
	x[i_x[k]] += 1.0/tweak;
      }

      multiscatter(n, m, config, instrument, surface, range, radius,
		   ext, ssa, x, ext_bscat_ratio, ext_air, ssa_air,
		   droplet_fraction, pristine_ice_fraction,
		   bscat_tweak, bscat_air_tweak);
      /* Loop through all the required measurements, which may be in
	 the multiple fields-of-view */
      for (l = 0; l < n_y; l++) {
	d_bscat_d_g[k][l] = (bscat_tweak[i_y[l]] - bscat[i_y[l]])
	  / (x[i_x[k]] - g[i_x[k]]);
	if (bscat_air) {
	  d_bscat_air_d_g[k][l] 
	    = (bscat_air_tweak[i_y[l]] - bscat_air[i_y[l]])
	    / (x[i_x[k]] - g[i_x[k]]);
	}
      }
      x[i_x[k]] = g[i_x[k]];
    }
  }


  /* JACOBIAN WITH RESPECT TO RADIUS */

  if (d_bscat_d_radius) {
    /* Copy radius into x */
    int i, k, l;
    for (i = 0; i < n; i++) {
      x[i] = radius[i];
    }

    /* Loop through all the x values to find the derivative for */
    for (k = 0; k < n_x; k++) {
      x[i_x[k]] *= exp_tweak;

      multiscatter(n, m, config, instrument, surface, range, x,
		   ext, ssa, g, ext_bscat_ratio, ext_air, ssa_air,
		   droplet_fraction, pristine_ice_fraction,
		   bscat_tweak, bscat_air_tweak);
      /* Loop through all the required measurements, which may be in
	 the multiple fields-of-view */
      for (l = 0; l < n_y; l++) {
	d_bscat_d_radius[k][l] = (bscat_tweak[i_y[l]] - bscat[i_y[l]])
	  / (x[i_x[k]] - radius[i_x[k]]);
	if (bscat_air) {
	  d_bscat_air_d_radius[k][l] 
	    = (bscat_air_tweak[i_y[l]] - bscat_air[i_y[l]])
	    / (x[i_x[k]] - radius[i_x[k]]);
	}
      }
      x[i_x[k]] = radius[i_x[k]];
    }
  }

  return MS_SUCCESS;
}


int
ms_jacobian_linear(
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
    const ms_real *pristine_ice_fraction,/* Fraction of ext due to pristine ice */
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
    )
{
  /* Only a numerical Jacobian is provided in the linear case */
  return ms_numerical_jacobian_linear(
           n, m, config, instrument, surface, 
	   range, radius, ext, ssa, g, ext_bscat_ratio,
	   ext_air, ssa_air,
	   droplet_fraction, pristine_ice_fraction,
	   n_x, n_y, i_x, i_y, 
	   bscat_out, bscat_air_out,
	   d_bscat_d_ext, d_bscat_d_ssa, d_bscat_d_g,
	   d_bscat_d_radius, d_bscat_d_ext_bscat_ratio,
	   d_bscat_air_d_ext, d_bscat_air_d_ssa, d_bscat_air_d_g,
	   d_bscat_air_d_radius);
}
