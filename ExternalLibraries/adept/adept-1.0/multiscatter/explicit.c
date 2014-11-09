/* explicit.c -- Algorithm treating scattering orders explicitly

   Copyright (C) 2007 Robin Hogan <r.j.hogan@reading.ac.uk> 

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ms.h"
#include "ms_switch.h"

static int ms_max_scattering_order = 4;

/* Function to be called recursively to calculate contributions from
   scattering up to order "norder" */
static
int
scatter(int n, int start_gate, int order, int norder, int output_dist,
	const ms_real *range,
	areal energy, areal width2, areal zeta2, areal cov,
	const areal *fs_fraction, const areal *theta2, 
	const ms_real *rho2range2, areal *bscat_ratio_out) {
  areal start_factor = 0.5;
  int i;
  
  for (i = start_gate; i < n; i++) {
    ms_real range_diff = fabs(range[i] - range[start_gate]);
    areal width2i = width2 + 2.0*cov*range_diff
      + zeta2*range_diff*range_diff;
    /* Add backscatter contribution */
    areal bscat_contrib = energy*start_factor
      *(1.0-exp(-rho2range2[i]/width2i));
    bscat_ratio_out[i] += bscat_contrib;
    /*
    if (order == 2) {
      ms_bscat_double[i] += bscat_contrib;
    }
    else {
      ms_bscat_multi[i] += bscat_contrib;
    }
    if (output_dist && i == n-1) {
      ms_increment_distribution(order, i, energy*start_factor, width2i);
    }
    */

    if (order < norder) {
      scatter(n, i, order+1, norder, output_dist, range,
	      energy*start_factor*fs_fraction[i], /* New energy */
	      width2i, /* New positional variance */
	      zeta2+theta2[i], /* New angular variance */
	      cov+zeta2*range_diff, /* New covariance */
	      fs_fraction, theta2, rho2range2, bscat_ratio_out);
    }
    start_factor = 1.0;
  }
  
  return MS_SUCCESS;
}

/* Use the photon variance-covariance approach to calculate the
   backscatter but treating each order of scattering explicitly. This
   is equivalent to Eloranta's formulation. */
int
ms_explicit(
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
    areal *bscat_air_out)
{
  int i;
  int output_dist = 0;
  ms_real rho_transmitter2 = instrument.rho_transmitter
    *instrument.rho_transmitter;
  areal transmittance = 1.0, optical_depth = 0.0;
  ms_real bscat_ext_ratio_air = 3.0 / (8.0 * MS_PI);
  ms_real rho_factor = 1.0/(1-exp(-instrument.rho_receiver[0]
				  *instrument.rho_receiver[0]/
				  (instrument.rho_transmitter
				  *instrument.rho_transmitter)));

  /* Initialise arrays; note that the ms_bscat* variables defined in
     multiscatter_qsa.c are also initialised */
  areal fs_fraction[n+1];
  areal theta2[n+1];
  ms_real rho2range2[n+1];
  ms_real drange[n];
  for (i = 0; i <= n; i++) {
    fs_fraction[i] = 0.0;
    theta2[i] = 0.0;
    rho2range2[i] = 0.0;
  }

  if (ms_init_intermediate_arrays(n) == MS_FAILURE) {
    return MS_FAILURE;
  }

  /* If we are to output the lateral photon distribution then
     initialise the relevant 2D arrays */
  if (config->options & MS_OUTPUT_DISTRIBUTION) {
    output_dist = 1;
    if (ms_init_distribution_arrays(n) == MS_FAILURE) {
      return MS_FAILURE;
    }
  }

  /* Pre-compute some convenient variables */
  for (i = 0; i < n; i++) {
    areal theta = instrument.wavelength / (MS_PI*radius[i]);
    drange[i] = ms_get_drange(i, n, range);
    theta2[i] = theta*theta;
    fs_fraction[i] = ext[i]*drange[i];
    rho2range2[i] = instrument.rho_receiver[0]*instrument.rho_receiver[0]
      * (range[i]-instrument.altitude)*(range[i]-instrument.altitude);
    bscat_out[i] = 0.0;
  }

  /* For double and higher-order scattering, call "scatter" which
     recursively calls itself for each subsequent order of scattering
     and each subsequent range gate */
  for (i = 0; i < n; i++) {
    ms_real range_abs = fabs(instrument.altitude-range[i]);
    scatter(n, i, 2, config->max_scattering_order, output_dist, 
	    range, fs_fraction[i],  /* Initial energy */
	    range_abs*range_abs*rho_transmitter2, /* Position variance */
	    rho_transmitter2+theta2[i], /* Angular variance */
	    rho_transmitter2*range_abs, /* Covariance */
	    fs_fraction, theta2, rho2range2,
	    bscat_out);
  }

  /* Work out the final backscatter */
  for (i = 0; i < n; i++) {
    areal od = ext[i]*drange[i]; /* optical depth of layer */
    if (od > 0.0) {
      transmittance = exp(-optical_depth)*(1-exp(-od))/od;
      optical_depth += od;
    }
    else {
      transmittance = exp(-optical_depth);
    }
    /* The single-scattering backscatter is the sum of the
       contribution from cloud/aerosol and air, multiplied by the
       transmittance squared. Note that we are assuming that if the
       air molecules can scatter then they will behave as Rayleigh
       scatterers. */
    if (ext_air) {
      areal transmittance_2way = transmittance*transmittance;
      /*
      ms_bscat_single[i] = transmittance_2way
	*(ext[i]/ext_bscat_ratio[i]
	  + ext_air[i]*bscat_ext_ratio_air);
      */
      if (bscat_air_out) {
	bscat_air_out[i] = (1.0+bscat_out[i]*rho_factor)
	  *transmittance_2way*ext_air[i]*bscat_ext_ratio_air;
	bscat_out[i] = (1.0+bscat_out[i]*rho_factor)
	  *transmittance_2way*ext[i]/ext_bscat_ratio[i];
      }
      else {
	bscat_out[i] = ms_bscat_single[i]*(1+bscat_out[i]*rho_factor);
      }
    }
    else {
      /*
      ms_bscat_single[i] = transmittance*transmittance*ext[i]
	/ext_bscat_ratio[i];
      */
      bscat_out[i] = ms_bscat_single[i]*(1+bscat_out[i]*rho_factor);
      if (bscat_air_out) {
	bscat_air_out[i] = 0.0;
      }
    }

    /* Scale the double and triple-and-higher component
       appropriately */
    /*
    ms_bscat_double[i] *= ms_bscat_single[i]*rho_factor;
    ms_bscat_multi[i] *= ms_bscat_single[i]*rho_factor;
    */
  }
  return MS_SUCCESS;
}

