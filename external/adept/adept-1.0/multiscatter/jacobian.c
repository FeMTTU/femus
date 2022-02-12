/* ms_jacobian.c -- Jacobian of multiscatter algorithm

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

#include "ms.h"

/* Calculate the Jacobian of the multiscatter algorithm by rerunning
   the algorithm for perturbations of the input values. This is slow
   but reliable. */
static
int
ms_numerical_jacobian(
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
    ms_real *droplet_fraction,/* Fraction of ext due to droplets, 0-1 */
    ms_real *ice_fraction,    /* Fraction of ext due to pristine ice */
    int n_x,  /* Number of state variables to compute the Jacobian of */
    int n_y,  /* Number of measurements to compute the Jacobian of */
    int *i_x, /* Index to the required state variables */
    int *i_y, /* Index to the required measurements */
    /* Output data: each Jacobian may be NULL if not required */
    ms_real *bscat,                  /* Measured backscatter, m-1 sr-1 */
    ms_real **d_ln_beta_d_ln_ext,    /* Jacobian with respect to ln(ext) */
    ms_real **d_ln_beta_d_ln_radius, /* Jacobian w.r.t. ln(radius) */
    ms_real *d_ln_beta_d_ln_ext_bscat_ratio /* ...w.r.t. ln(ext_bscat_ratio) */
    )
{
  ms_real tweak = 1000;
  ms_real exp_tweak = exp(1.0/tweak);
  ms_real bscat_tweak[n];
  ms_real x[n];
  int i;
  /* Ordinary run */
  multiscatter(n, m, config, instrument, surface, range, radius,
	       ext, ssa, g, ext_bscat_ratio, ext_air, ssa_air,
	       droplet_fraction, ice_fraction,
	       bscat, NULL);

  /* Set each element of each array to zero. */
  for (i = 0; i < n; i++) {
    bscat_tweak[i] = 0.0;
    x[i] = 1.0;
  }


  if (d_ln_beta_d_ln_ext_bscat_ratio) {
    int l;
    for (l = 0; l < n_y; l++) {
      if (ext[i_y[l]] > 0.0) {
	d_ln_beta_d_ln_ext_bscat_ratio[l] = -(ms_bscat_single[i_y[l]]
		      + ms_bscat_double[i_y[l]] + ms_bscat_multi[i_y[l]])
	  / bscat[i_y[l]];
      }
      else {
	d_ln_beta_d_ln_ext_bscat_ratio[l] = 0.0;
      }
    }
    /*
    int k;
    for (k = 0; k < n_x; k++) {
      d_ln_beta_d_ln_ext_bscat_ratio[k] = -(ms_bscat_single[i_x[k]]
		 + ms_bscat_double[i_x[k]] + ms_bscat_multi[i_x[k]])/bscat[i_x[k]];
    }
    */
  }
  if (d_ln_beta_d_ln_ext) {
    /* Copy ext into x */
    int i, k, l;
    for (i = 0; i < n; i++) {
      x[i] = ext[i];
    }
    /* Loop through all the x values to find the derivative for */
    for (k = 0; k < n_x; k++) {
      x[i_x[k]] *= exp_tweak;

      multiscatter(n, m, config, instrument, surface, range, radius,
		   x, ssa, g, ext_bscat_ratio, ext_air, ssa_air,
		   droplet_fraction, ice_fraction,
		   bscat_tweak, NULL);
      for (l = 0; l < n_y; l++) {
	d_ln_beta_d_ln_ext[l][k] = tweak * (bscat_tweak[i_y[l]]/bscat[i_y[l]]
					    - 1.0);
      }
      x[i_x[k]] = ext[i_x[k]];
    }
  }
  if (d_ln_beta_d_ln_radius) {
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
		   droplet_fraction, ice_fraction,
		   bscat_tweak, NULL);
      for (l = 0; l < n_y; l++) {
	d_ln_beta_d_ln_radius[l][k] = tweak * (bscat_tweak[i_y[l]]/bscat[i_y[l]] - 1.0);
      }
      x[i_x[k]] = radius[i_x[k]];
    }
  }
  return MS_SUCCESS;
}

static
int
ms_analytic_jacobian(
    /* Input data */
    int n,                    /* Number of range gates in profile */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    ms_real *range,           /* Height of each range gate, metres */
    ms_real *radius,          /* Cloud/aerosol equivalent radius, microns */
    ms_real *ext,             /* Cloud/aerosol extinction coefficient, m-1 */
    ms_real *ext_bscat_ratio, /* Cloud/aerosol ext./backscatter ratio, sr */
    ms_real *ext_air,         /* Air ext. coefft, m-1 (NULL for vacuum) */
    ms_real *droplet_fraction,/* Fraction of ext due to droplets, 0-1 */
    ms_real *ice_fraction,    /* Fraction of ext due to pristine ice */
    int n_x,  /* Number of state variables to compute the Jacobian of */
    int n_y,  /* Number of measurements to compute the Jacobian of */
    int *i_x, /* Index to the required state variables */
    int *i_y, /* Index to the required measurements */
    /* Output data */
    ms_real *bscat,               /* Measured backscatter, m-1 sr-1 */
    ms_real **d_ln_beta_d_ln_ext,    /* Jacobian with respect to ln(ext) */
    ms_real **d_ln_beta_d_ln_radius, /* Jacobian w.r.t. ln(radius) */
    ms_real *d_ln_beta_d_ln_ext_bscat_ratio /* ...w.r.t. ln(ext_bscat_ratio) */
    )
{
  /* Variance of transmitter */
  ms_real mu2_transmitter = instrument.rho_transmitter
    *instrument.rho_transmitter; /* radians^2 */

  /* Variance of receiver: this is a top-hat distribution */
  ms_real mu2_receiver = 0.5 * instrument.rho_receiver[0]
    *instrument.rho_receiver[0];  /* radians^2 */
  ms_real fov_factor = 1.0 / (1.0 - exp(-instrument.rho_receiver[0]
					*instrument.rho_receiver[0]
				     / (instrument.rho_transmitter
					*instrument.rho_transmitter)));
  int range_sign_factor = range[2] < range[1] ? -1 : 1;

  ms_real self_extinction_factor = 0.5;
  ms_real self_scatter_factor = 0.5;

  ms_real drange = fabs(range[2]-range[1]);

  if (config->options & MS_CRUDE_INTEGRATION) {
    self_scatter_factor = 1.0;
  }
  else if (config->options & MS_NO_MULTISCAT_WITHIN_GATE) {
    self_scatter_factor = 0.0;
  }

  if (config->options & MS_CRUDE_OPTICAL_DEPTH) {
    self_extinction_factor = 1.0;
  }

  ms_small_angle(n, config, instrument, surface, range, radius,
		 ext, ext_bscat_ratio, ext_air,
		 droplet_fraction, ice_fraction,
		 bscat, NULL);

  if (d_ln_beta_d_ln_ext_bscat_ratio) {
    int l;
    for (l = 0; l < n_y; l++) {
      if (ext[i_y[l]] > 0.0) {
	d_ln_beta_d_ln_ext_bscat_ratio[l] = -(ms_bscat_single[i_y[l]]
		      + ms_bscat_double[i_y[l]] + ms_bscat_multi[i_y[l]])
	  / bscat[i_y[l]];
      }
      else {
	d_ln_beta_d_ln_ext_bscat_ratio[l] = 0.0;
      }
    }
  }
  if (d_ln_beta_d_ln_ext && d_ln_beta_d_ln_radius) {
    int k, l;
    /* Loop through all the x values to find the derivative for */
    for (k = 0; k < n_x; k++) {
      ms_real theta = instrument.wavelength / (MS_PI*radius[i_x[k]]);
      ms_real theta2 = theta*theta;
      for (l = 0; l < n_y; l++) {
	  ms_real total_range = range_sign_factor 
	    * (range[i_y[l]]-instrument.altitude);
	  ms_real total_range2 = total_range * total_range;
	  ms_real delta_range = range_sign_factor 
	    * (range[i_y[l]]-range[i_x[k]]);
	  ms_real delta_range2 = delta_range * delta_range;
	/* Single scattering */
	d_ln_beta_d_ln_radius[l][k] = 0.0;

	if (i_x[k] > i_y[l]) {
	  /* Observation pixel closer to lidar than state pixel: no
	     effect of state on observation and Jacobian = 0 */
	  d_ln_beta_d_ln_ext[l][k] = 0.0;
	}
	else if (i_x[k] < i_y[l]) {
	  ms_real denominator = theta2*delta_range2/total_range2+mu2_transmitter;
	  ms_real exp_term = exp(-2.0*mu2_receiver / denominator);
	  /* Observation pixel is further from lidar than state pixel:
	     effect due to extinction and multiple scattering */
	  d_ln_beta_d_ln_ext[l][k] =
	    /* extinction: */
	    -2.0 * ext[i_x[k]] * drange
	    /* double scattering: */
	    + (ms_bscat_single[i_y[l]] + ms_bscat_air_single[i_y[l]])
	    * fov_factor * drange * ext[i_x[k]]
	    * (1 - exp_term) / bscat[i_y[l]];
	  /* Radius term due to double scattering */
	  d_ln_beta_d_ln_radius[l][k] = 
	    (ms_bscat_single[i_y[l]] + ms_bscat_air_single[i_y[l]])
	    * fov_factor * drange * ext[i_x[k]] / bscat[i_y[l]]
	    * exp_term
	    * 2.0*mu2_receiver/(denominator*denominator)
	    * 2.0*theta2*delta_range2/total_range2;
	}
	else {
	  /* Observation and state pixel are the same: effect due to
	     extinction, multiple scattering and the ext/bscat
	     ratio */ 
	  d_ln_beta_d_ln_ext[l][k] = (
	    /* Single scattering: */
	    ms_bscat_single[i_x[k]] + ms_bscat_double[i_x[k]]
	                            + ms_bscat_multi[i_x[k]]
	    -self_extinction_factor * 2.0 * ext[i_x[k]] 
	                            * bscat[i_y[l]] * drange
	    /* Double scattering */
	    +self_scatter_factor * (ms_bscat_single[i_y[l]]
				    + ms_bscat_air_single[i_y[l]])
	    * fov_factor * fov_factor * drange * ext[i_x[k]]) / bscat[i_y[l]];
	}
      }
    }
  }
  return MS_SUCCESS;
}

int
ms_jacobian(
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
    ms_real *droplet_fraction,/* Fraction of ext due to droplets, 0-1 */
    ms_real *ice_fraction,    /* Fraction of ext due to pristine ice */
    int n_x,  /* Number of state variables to compute the Jacobian of */
    int n_y,  /* Number of measurements to compute the Jacobian of */
    int *i_x, /* Index to the required state variables */
    int *i_y, /* Index to the required measurements */
    /* Output data */
    ms_real *bscat_out,           /* Measured backscatter, m-1 sr-1 */
    ms_real **d_ln_beta_d_ln_ext,    /* Jacobian with respect to ln(ext) */
    ms_real **d_ln_beta_d_ln_radius, /* Jacobian w.r.t. ln(radius) */
    ms_real *d_ln_beta_d_ln_ext_bscat_ratio /* ...w.r.t. ln(ext_bscat_ratio) */
    )
{
  int status;

  if (config->options & MS_NUMERICAL_JACOBIAN) {
    status = ms_numerical_jacobian(
           n, m, config, instrument, surface, 
	   range, radius, ext, ssa, g, ext_bscat_ratio,
	   ext_air, ssa_air,
	   droplet_fraction, ice_fraction,
	   n_x, n_y, i_x, i_y, 
	   bscat_out,
	   d_ln_beta_d_ln_ext, d_ln_beta_d_ln_radius, 
	   d_ln_beta_d_ln_ext_bscat_ratio);
  }
  else {
    if (config->wide_angle_algorithm != MS_WIDE_ANGLE_NONE) {
      fprintf(stderr, "Error: analytic Jacobian for time-dependent two-stream method is not yet available\n");
      return MS_FAILURE;
    }
    status = ms_analytic_jacobian(
           n, config, instrument, surface,
	   range, radius, ext, ext_bscat_ratio, ext_air,
	   droplet_fraction, ice_fraction,
	   n_x, n_y, i_x, i_y, 
	   bscat_out,
	   d_ln_beta_d_ln_ext, d_ln_beta_d_ln_radius, 
	   d_ln_beta_d_ln_ext_bscat_ratio);
  }
  return status;
}
