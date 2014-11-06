/* multiscatter.c -- Interface to multiple scattering algorithms

   Copyright (C) 2006-2010 Robin Hogan <r.j.hogan@reading.ac.uk> 

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


   The algorithm is implemented in ANSI C, but a Fortran interface is
   provided.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* For details of how to call the multiple scattering algorithms, see
   the comments in multiscatter.h, called from ms.h: */
#include "ms.h"
#include "ms_switch.h"

/* Perform the multiple scattering calculation, using which ever
   combination of algorithms is appropriate depending on multiple
   scattering regime */
int
multiscatter(
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
    areal *bscat_air_out)   /* Measured backscatter of air, m-1 sr-1 */
{
  static char have_warned_about_calibration_factor = 0;
  ms_instrument instrument_1fov = instrument;
  ms_real ss_multiplier = 1.0;
  int status = MS_SUCCESS;
  int i;
  int ifov;
  int require_bscat_air = (bscat_air_out != NULL);

  /* This variable contains the same instrument settings, except for
     having only one field-of-view; it is modified when we loop over
     the fields-of-view */
  instrument_1fov.nfov = 1;

  /* Work out single-scattering factor: this ensures that the first
     field-of-view is properly calibrated and that the backscatter in
     the subsequent fields-of-view are calibrated to the first such
     that the backscatter is proportional to the energy received. */
  if (instrument.receiver_type == MS_TOP_HAT) {
    ss_multiplier
      = 1.0/(1.0-exp(-instrument.rho_receiver[0]*instrument.rho_receiver[0]/
		     (instrument.rho_transmitter*instrument.rho_transmitter)));
  }
  else {
    ss_multiplier
      = 1.0 + instrument.rho_transmitter*instrument.rho_transmitter
      /(instrument.rho_receiver[0]*instrument.rho_receiver[0]);
  }
  config->ss_multiplier = ss_multiplier;

  if (!(config->options & MS_QUIET) && ss_multiplier != 1.0
      && !have_warned_about_calibration_factor) {
    fprintf(stderr, "Factor to ensure calibration of first FOV: %g\n",
	    ss_multiplier);
    have_warned_about_calibration_factor = 1;
  }


  /* FIRST DO SMALL-ANGLE SCATTERING CALCULATION */

  if (config->small_angle_algorithm != MS_SINGLE_AND_SMALL_ANGLE_NONE) {
    /* For small-angle scattering we loop over the fields-of-view */
    for (ifov = 0; ifov < instrument.nfov; ifov++) {
      /* Set the receiver field of view appropriately */
      ms_real rho_receiver = instrument.rho_receiver[ifov];
      instrument_1fov.rho_receiver = &rho_receiver;
      
      if (config->small_angle_algorithm == MS_SINGLE_SCATTERING) {
	/* If no forward lobe (radar case) then we have no small-angle
	   multiple scattering contribution and we simply calculate the
	   single scattering here, with wide-angle scattering later. */
	status = ms_singlescatter(n, instrument_1fov, surface,
				  range, ext, ext_bscat_ratio,
				  ext_air,
				  bscat_out + ifov*m,
				  bscat_air_out + ifov*m*require_bscat_air);
      }
      else if (config->small_angle_algorithm != MS_SMALL_ANGLE_PVC_EXPLICIT) {
	/* Perform small-angle multiple-scattering calculation using an
	   efficient algorithm */
	if (config->small_angle_algorithm == MS_SMALL_ANGLE_PVC_FAST_LAG) {
	  fprintf(stderr, "Error: Code computing lag is not included in this test version\n");
	  return MS_INPUT_ERROR;
	}
	else if (config->small_angle_algorithm == MS_SMALL_ANGLE_PVC_FAST) {
	  status = ms_fast_small_angle(n, 
		       instrument_1fov, surface, 
		       range, radius, ext, ext_bscat_ratio, ext_air,
		       droplet_fraction, pristine_ice_fraction,
		       bscat_out + ifov*m,
		       bscat_air_out + ifov*m*require_bscat_air);
	}
	else if (config->small_angle_algorithm == MS_SMALL_ANGLE_PVC_ORIGINAL) {
	  fprintf(stderr, "Error: Original small-angle code not included in this test version\n");
	  return MS_INPUT_ERROR;
	}
	else {
	  fprintf(stderr, "Error: no small-angle algorithm with code %d\n",
		  config->small_angle_algorithm);
	  return MS_INPUT_ERROR;
	}
      }
      else {
	/* Perform small-angle multiple-scattering calculation using an
	   explicit (slow) algorithm */
	status = ms_explicit(n, 
		     config, instrument_1fov, surface, 
		     range, radius, ext, ext_bscat_ratio, ext_air,
		     bscat_out + ifov*m,
		     bscat_air_out + ifov*m*require_bscat_air);
      }
      
      /* To ensure the correct calibration, we may need to scale the
	 wider fields of view - here we multiply by the multiplier for
	 the first field-of-view, and divide by the multiplier
	 calculated for the current field-of-view.  At the wide-angle
	 scattering stage, we simply multiply by the multiplier for
	 the first field-of-view.  This ensures that for a given
	 field-of-view, the narrow-angle and wide-angle contributions
	 are correctly scaled relative to each other (normally one
	 would only multiply the wide-angle component by the
	 multiplier for the current FOV), and relative to the first
	 field-of-view */
      if (ifov > 0) {
	ms_real fov_ss_multiplier;
	areal* bscat_cur = bscat_out+ifov*m;
	areal* bscat_air_cur = bscat_air_out+require_bscat_air*ifov*m;
	if (instrument.receiver_type == MS_TOP_HAT) {
	  fov_ss_multiplier
	    = ss_multiplier*(1.0-exp(-instrument.rho_receiver[ifov]
				     *instrument.rho_receiver[ifov]/
				     (instrument.rho_transmitter
				      *instrument.rho_transmitter)));
	}
	else {
	  fov_ss_multiplier = ss_multiplier
	    / (1.0 + instrument.rho_transmitter*instrument.rho_transmitter
	       /(instrument.rho_receiver[ifov]*instrument.rho_receiver[ifov]));
	}
	for (i = 0; i < n; i++, bscat_cur++) {
	  (*bscat_cur) *= fov_ss_multiplier;
	}
	if (bscat_air_cur) {
	  for (i = 0; i < n; i++, bscat_air_cur++) {
	    (*bscat_air_cur) *= fov_ss_multiplier;
	  }
	}
	if (!(config->options & MS_QUIET) && ss_multiplier != 1.0) {
	  fprintf(stderr, "   Equivalent factor for FOV %d: %g\n",
		  ifov+1, ss_multiplier/fov_ss_multiplier);
	}
      }
    } /* Loop over fields-of-view */

    if (config->options & MS_ANNULAR_DETECTORS) {
      /* Need to subtract backscatter contributions from each other */
      for (ifov = instrument.nfov-1; ifov > 0; ifov--) {
	int i;
	int i_inner = (ifov-1)*m;
	int i_outer = ifov*m;
	for (i = 0; i < m; i++, i_inner++, i_outer++) {
	  bscat_out[i_outer] -= bscat_out[i_inner];
	}
	if (bscat_air_out) {
	  i_inner = (ifov-1)*m;
	  i_outer = ifov*m;
	  for (i = 0; i < m; i++, i_inner++, i_outer++) {
	    bscat_air_out[i_outer] -= bscat_air_out[i_inner];
	  }
	}
      }
    }
  } /* If perform small-angle calculation */

  /* Return if status is non-zero (indicating that an error
     occurred) or if ssa or g are not set (indicating that the
     wide-angle calculation is not to be performed). */
  if (status != MS_SUCCESS || !ssa || !g
      || config->wide_angle_algorithm == MS_WIDE_ANGLE_NONE) {
    return status;
  }

  if (config->small_angle_algorithm == MS_SINGLE_AND_SMALL_ANGLE_NONE) {
    /* Reset backscatter to zero, since wide-angle calculation
       increments existing values */
    int i;
    for (i = 0; i < m*instrument.nfov; i++) {
      bscat_out[i] = 0.0;
    }
  }

  /* NOW DO WIDE-ANGLE SCATTERING CALCULATION */
  /*
  return ms_wide_angle_regrid(n, m,
			      config, instrument, surface, range, radius,
			      ext, ssa, g, ext_air, ssa_air,
			      bscat_out);
  */
  return ms_wide_angle(n, m,
			      config, instrument, surface, range, radius,
			      ext, ssa, g, ext_air, ssa_air,
			      bscat_out);
  
}

