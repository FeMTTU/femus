/* small_angle.c -- Small-angle multiple scattering using Hogan (2006) algorithm

   Copyright (C) 2004-2007 Robin Hogan <r.j.hogan@reading.ac.uk> 

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


   This file contains an algorithm for efficient calculation of the
   lidar backscatter profile in the presence of multiple scattering.
   It uses Eloranta's formulation for double scattering and
   parameterises the photon distribution for estimating scattering at
   higher orders in such a way that O(N^2) efficiency is ensured,
   where N is the number of range gates.

   There is a certain resolution dependence to this algorithm so it is
   best to choose the resolution such that the optical depth within
   any given layer is always considerably less than one.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ms.h"

/* The following intermediate variables are used to describe the
   photon distribution at each lidar range gate. They are available to
   external functions after the algorithm has been run. NOTE: they are
   calculated using the "equivalent-medium theorem", so use doubled
   extinction/scattering coefficients on the outward journey. */

/* Energy in outgoing photons that have been forward-scattered once,
   relative to the energy in the unscattered beam: */
ms_real *ms_E_once = NULL;

/* Second moment of the distribution of photons that have been
   forward-scattered once, expressed in terms of the angular deviation
   from the lidar axis. */
ms_real *ms_Emu2_once = NULL; /* radians^2 */

/* Energy in outgoing photons that have been forward-scattered more
   than once, relative to the energy in the unscattered beam:*/ 
ms_real *ms_E_multi = NULL;

/* Second moment of the distribution of photons that have been
   forward-scattered more than once, expressed in terms of the angular
   deviation from the lidar axis: */
ms_real *ms_Emu2_multi = NULL;  /* radians^2 */

/* Second moment of distribution of photon pointing angles relative to
   direction to lidar, for the once scattered and multiply scattered
   photons: */
ms_real *ms_Pmu2_once = NULL;   /* radians^2 */
ms_real *ms_Pmu2_multi = NULL;  /* radians^2 */

/* Covariance of the photon pointing angles and their angular
   distribution relative to the lidar: */
ms_real *ms_EPcov_once = NULL;  /* radians^2 */
ms_real *ms_EPcov_multi = NULL; /* radians^2 */

/* Single-scattering only backscatter, m-1 sr-1 */
ms_real *ms_bscat_single = NULL;

/* Double-scattering only backscatter, m-1 sr-1 */
ms_real *ms_bscat_double = NULL;

/* Triple and higher scattering, m-1 sr-1 */
ms_real *ms_bscat_multi = NULL;

/* Single-scattering only backscatter, m-1 sr-1 */
ms_real *ms_bscat_air_single = NULL;

/* Double-scattering only backscatter, m-1 sr-1 */
ms_real *ms_bscat_air_double = NULL;

/* Triple and higher scattering, m-1 sr-1 */
ms_real *ms_bscat_air_multi = NULL;

/* Bit-field containing algorithm options */
int ms_options = 0;

/* If the forward scattering angle is above this value (in radians)
   then the radiation is deemed to have escaped and cannot be returned
   to the receiver. This is to prevent problems associated with
   aerosols. Note that double scattering, which is unaffected by this
   problem, is still calculated. */
ms_real ms_max_theta = 0.1;

/* Number of gates currently allocated in the intermediate
   variables: */
static int ms_num_gates_allocated = 0;

/* Number of gates defined in the last call to the algorithm: */
static int ms_num_gates_defined = 0;


/* Intermediate arrays are allocated dynamically and can be read by
   external functions after the algorithm has been run. For speed,
   subsequent runs of the algorithm will use the same memory but new
   memory will be allocated if there is not enough space for the new
   profile. Whether or not new memory is allocated, the arrays are
   always set to zero. Note that this function is called automatically
   by multiscater(). It returns MS_FAILURE if there is a problem
   allocating the memory, MS_SUCCESS otherwise. */
int
ms_init_intermediate_arrays(int num_gates)
{
  /* Macro for allocating or reallocating memory for a single array */
#define ALLOCATE(target, num_bytes) \
  if ((target = (ms_real*)realloc(target, num_bytes)) == NULL) {	\
      target = NULL; return MS_FAILURE; }

  int i;
  if (ms_num_gates_allocated < num_gates) {
    /* Number of gates allocated is less than number requested (which
       may be zero if this function has not been called): (re)allocate
       memory. */
    int num_bytes = num_gates * sizeof(ms_real);
    ms_num_gates_allocated = 0; /* In case error occurs */
    ALLOCATE(ms_E_once, num_bytes);
    ALLOCATE(ms_Emu2_once, num_bytes);
    ALLOCATE(ms_E_multi, num_bytes);
    ALLOCATE(ms_Emu2_multi, num_bytes);
    ALLOCATE(ms_Pmu2_once, num_bytes);
    ALLOCATE(ms_Pmu2_multi, num_bytes);
    ALLOCATE(ms_EPcov_once, num_bytes);
    ALLOCATE(ms_EPcov_multi, num_bytes);
    ALLOCATE(ms_bscat_single, num_bytes);
    ALLOCATE(ms_bscat_double, num_bytes);
    ALLOCATE(ms_bscat_multi, num_bytes);
    ALLOCATE(ms_bscat_air_single, num_bytes);
    ALLOCATE(ms_bscat_air_double, num_bytes);
    ALLOCATE(ms_bscat_air_multi, num_bytes);
  }
  ms_num_gates_allocated = num_gates;
  /* Set each element of each array to zero. */
  for (i = 0; i < num_gates; i++) {
    ms_E_once[i] = 0.0;
    ms_Emu2_once[i] = 0.0;
    ms_E_multi[i] = 0.0;
    ms_Emu2_multi[i] = 0.0;
    ms_Pmu2_once[i] = 0.0;
    ms_Pmu2_multi[i] = 0.0;
    ms_EPcov_once[i] = 0.0;
    ms_EPcov_multi[i] = 0.0;
    ms_bscat_single[i] = 0.0;
    ms_bscat_double[i] = 0.0;
    ms_bscat_multi[i] = 0.0;
    ms_bscat_air_single[i] = 0.0;
    ms_bscat_air_double[i] = 0.0;
    ms_bscat_air_multi[i] = 0.0;
  }
  ms_num_gates_defined = num_gates;
  return MS_SUCCESS;
}

/* Optional: deallocate memory used by intermediate arrays. */
void
ms_free_intermediate_arrays()
{
  if (ms_num_gates_allocated > 0) {
    free(ms_E_once);
    free(ms_Emu2_once);
    free(ms_E_multi);
    free(ms_Emu2_multi);
    free(ms_Pmu2_once);
    free(ms_Pmu2_multi);
    free(ms_EPcov_once);
    free(ms_EPcov_multi);
    free(ms_bscat_single);
    free(ms_bscat_double);
    free(ms_bscat_multi);
    free(ms_bscat_air_single);
    free(ms_bscat_air_double);
    free(ms_bscat_air_multi);
  }
  ms_num_gates_allocated = 0;
  ms_num_gates_defined = 0;
}


/* C interface: this function runs the algorithm. The data are entered
   as pointers to arrays of length n, and the function writes to the
   memory pointed to by bscat_out. See the header file for the
   meanings and units of the input arguments. The function returns
   MS_FAILURE if memory cannot be allocated for the intermediate
   arrays. */
int
ms_small_angle(
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
    ms_real *bscat_air_out)   /* Measured backscatter of air, m-1 sr-1 */
{
  /* Variance of transmitter */
  ms_real mu2_transmitter
    = instrument.rho_transmitter*instrument.rho_transmitter; /* radians^2 */

  /* Variance of receiver: this is a top-hat distribution */
  ms_real mu2_receiver = 0.5
    * instrument.rho_receiver[0]*instrument.rho_receiver[0];  /* radians^2 */
  ms_real fov_factor = 1.0
    / (1.0 - exp(-instrument.rho_receiver[0]*instrument.rho_receiver[0]
		 / (instrument.rho_transmitter*instrument.rho_transmitter)));

  /* Optical depth including air */
  ms_real total_optical_depth = 0.0;
  ms_real transmission_factor = 1.0;

  /* Optical depth of cloud/aerosol only, for forward-scattering
     calculations */
  ms_real fs_optical_depth = 0.0;

  ms_real last_ext = 0.0, last_ext_air = 0.0;

  ms_real bscat_ext_ratio_air = 3.0 / (8.0 * MS_PI);

  int range_sign_factor = range[2] < range[1] ? -1 : 1;

  int i, j;


  /* Allocate and set-to-zero the intermediate arrays */
  if (ms_init_intermediate_arrays(n) == MS_FAILURE) {
    return MS_FAILURE;
  }

  /* Assume that for wavelengths longer than 1 micron the gaseous
     extinction is entirely due to absorption rather than Rayleigh
     scattering */
  if (instrument.wavelength > 1e-6) {
    bscat_ext_ratio_air = 0.0;
  }

  /* Treat each lidar range gate in turn */
  for (i = 0; i < n; i++) {
    ms_real drange = ms_get_drange(i, n, range);

    ms_real bscat_air = ext_air[i] * bscat_ext_ratio_air;
    ms_real trapez_int_factor = 0.5;

    if (config->options & MS_NO_MOLECULAR_EXTINCTION) {
      if (config->options & MS_SIMPLE_OPTICAL_DEPTH) {
	total_optical_depth += 0.5 * (ext[i]+last_ext)*drange;
      }
      else if (config->options & MS_CRUDE_OPTICAL_DEPTH) {
	total_optical_depth += ext[i]*drange;
      }
      else {
	total_optical_depth += last_ext*drange;
	if (ext[i] > 0.0) {
	  transmission_factor = (1-exp(-2.0*ext[i]*drange))
	    / (2.0*ext[i]*drange);
	}
	else {
	  transmission_factor = 1.0;
	}
      }
    }
    else {
      if (config->options & MS_SIMPLE_OPTICAL_DEPTH) {
	total_optical_depth += 0.5 * (ext[i]+last_ext
				      +ext_air[i]+last_ext_air)*drange;
      }
      else if (config->options & MS_CRUDE_OPTICAL_DEPTH) {
	total_optical_depth += (ext[i]+ext_air[i])*drange;
      }
      else {
	total_optical_depth += (last_ext+last_ext_air)*drange;
	if (ext[i] > 0.0) {
	  transmission_factor = (1-exp(-2.0*(ext[i]+ext_air[i])*drange))
	    / (2.0*(ext[i]+ext_air[i])*drange);
	}
	else {
	  transmission_factor = 1.0;
	}
      }
    }

    if (config->options & MS_CRUDE_INTEGRATION) {
      fs_optical_depth += ext[i]*drange;
      trapez_int_factor = 1.0;
    }
    else if (config->options & MS_NO_MULTISCAT_WITHIN_GATE) {
      fs_optical_depth += last_ext*drange;
      trapez_int_factor = 0.0;
    }
    else {
      fs_optical_depth += 0.5 * (ext[i]+last_ext)*drange;
    }

    last_ext = ext[i];
    last_ext_air = ext_air[i];

    /* Multiple-scattering calculation only carried out for this gate
       if there is a non-zero cloud/aerosol extinction (molecular
       scattering pattern does not lead to significant multiple
       scattering) */
    if (ext[i] > 0.0) {
      /* Assume that half the extinguished radiation is in the forward
	 lobe, but use the Katsev et al. theorem and use an
	 equivalent medium where the extinction is twice reality on
	 the outward journey but zero on the return journey */
      ms_real fraction_forward_scattered = ext[i]*drange;
      /* Width of forward scattered lobe */
      ms_real theta = instrument.wavelength / (MS_PI*radius[i]);
      ms_real theta2 = theta*theta;

      /* Crude way to deal with thin aerosol layers with small
	 particles that lead to very non-Gaussian distributions when
	 there is also a cloud in the profile. */
      if (config->options & MS_WIDE_ANGLE_CUTOFF && theta > ms_max_theta) {
	fraction_forward_scattered = 0.0;
      }

      /* We now work out the new contributions to the once- and
         multi-scattered distributions.  The only source for the
         once-scattered distribution is the unscattered transmitter
         beam. The multi-scattered distribution has sources from the
         once-scattered distribution as well as rescattered photons
         from the multi-scattered distribution. */

      /* Calculate the contribution of the forward scattered energy at
	 gate i to all gates further away. */
      for (j = i; j < n; j++) {
	ms_real total_range = range_sign_factor * (range[j]-instrument.altitude);
	ms_real total_range2 = total_range * total_range;
	ms_real delta_range = range_sign_factor * (range[j]-range[i]);
	ms_real delta_range2 = delta_range * delta_range;
	ms_real weight = trapez_int_factor * fraction_forward_scattered;
	if (config->options & MS_APPROXIMATE_EXPONENTIAL) {
	  ms_bscat_double[j] += trapez_int_factor * ext[i]
	    * mu2_receiver/(theta2*delta_range2/total_range2 + mu2_transmitter);
	}
	else if (!(config->options & MS_CRUDE_DOUBLE_SCATTERING)) {
	  /* Double scattering contribution using Eloranta's method;
	     note that ms_bscat_double is rescaled later. */
	  /* This call to exp() slows the whole algorithm down */
	  ms_bscat_double[j] += trapez_int_factor * ext[i]
	    * exp(-2.0*mu2_receiver
		  / (theta2*delta_range2/total_range2 + mu2_transmitter));
	}
	/* Weighted average of the moments of the once-scattered
	   distribution */
	ms_E_once[j] += weight;
	ms_Emu2_once[j] += weight*(mu2_transmitter*total_range2 + theta2*delta_range2);
	ms_Pmu2_once[j] += weight*(mu2_transmitter+theta2);
	ms_EPcov_once[j] += weight
	  * (mu2_transmitter*total_range + delta_range*theta2);

	if (ms_E_once[i] > 0.0) {
	  /* Weighted average of the moments of the multiply-scattered
	     distribution */
	  ms_Emu2_multi[j] += weight 
	    * ( (ms_Emu2_once[i] + ms_Emu2_multi[i])
		+ (ms_Pmu2_once[i] + ms_Pmu2_multi[i]
		   + theta2*(ms_E_once[i]+ms_E_multi[i])) * delta_range2
		+ 2.0*(ms_EPcov_once[i]+ms_EPcov_multi[i])* delta_range);
	  
	  ms_Pmu2_multi[j] += weight
	    * (ms_Pmu2_once[i] + ms_Pmu2_multi[i]
	       + theta2*(ms_E_once[i]+ms_E_multi[i]));

	  ms_EPcov_multi[j] += weight
	    * (ms_EPcov_once[i] + ms_EPcov_multi[i]
	       + (theta2*(ms_E_once[i]+ms_E_multi[i])
		  + ms_Pmu2_once[i] + ms_Pmu2_multi[i]) * delta_range);

	  ms_E_multi[j] += weight*(ms_E_once[i] + ms_E_multi[i]);
	}
	trapez_int_factor = 1.0;
      }
    }

    if ((ext_air && ext_air[i] > 0.0) || ext[i] > 0.0) {
      ms_real integral_ratio = 0.0;
      ms_real two_way_transmission = transmission_factor
	* exp(-2.0 * total_optical_depth);
      ms_real total_range2 = (range[i]-instrument.altitude)
	* (range[i]-instrument.altitude);

      ms_real afactor = 1.0;

      if (ext[i] > 0.0 && droplet_fraction && pristine_ice_fraction 
	  && (droplet_fraction[i] > 0.0 || pristine_ice_fraction[i] > 0.0)) {
	/* The particle phase functions in the near-180-degree
	   direction are not isotropic, so calculate the factor by
	   which the backscatter should be multiplied to account for
	   this effect */
	ms_real Escale = 1.0/(ms_E_once[i]+ms_E_multi[i]);
	ms_real Range = fabs(range[i]-instrument.altitude);
	afactor = ms_anisotropic_factor(
	     /* Properties of forward scattered distribution */
	     (ms_Emu2_once[i]+ms_Emu2_multi[i])*Escale, /* Spatial variance */
	     (ms_Pmu2_once[i]+ms_Pmu2_multi[i])*Escale, /* Angular variance */
	     (ms_EPcov_once[i]+ms_EPcov_multi[i])*Escale, /* Covariance */
	     instrument.wavelength*instrument.wavelength
	     / (MS_PI*MS_PI*radius[i]*radius[i]), /* Theta squared */
	     Range,  Range*Range*instrument.rho_receiver[0]
	     *instrument.rho_receiver[0], /* Receiver FOV squared */
	     ext[i]/ext_bscat_ratio[i], bscat_air, /* Backscatter weights */
	     droplet_fraction[i], pristine_ice_fraction[i]);
      }

      /* Calculate the backscatter due to single scattering */
      ms_bscat_single[i] = two_way_transmission * ext[i] / ext_bscat_ratio[i];
      ms_bscat_air_single[i] = two_way_transmission * bscat_air;

      /* Integrate to calculate backscatter due to higher-order
	 scatterings, divided by the single scatter value */
#define INT_GAUSSIAN(E, K, Emu2, Kmu2) ((E)*(K)/((Emu2)/(E)+(Kmu2)/(K)))
#define INT_TOPHAT(E, K, Emu2, Kmu2) ((E)*(K)) \
      *(1.0-exp(-2.0*(Kmu2)*(E)/((K)*(Emu2))))/(2.0*(Kmu2)/(K))

      /* Calculate the backscatter due to double scattering */
      if (config->options & MS_APPROXIMATE_EXPONENTIAL) {
	ms_bscat_air_double[i] = ms_bscat_air_single[i] * fov_factor 
	  * ms_bscat_double[i]*drange;
	/* Now set ms_bscat_double to its proper value */
	ms_bscat_double[i] = ms_bscat_single[i] * fov_factor 
	  * ms_bscat_double[i]*drange;
      }
      if (!(config->options & MS_CRUDE_DOUBLE_SCATTERING)) {
	/* Use Eloranta's double scattering formulation */
	/* Note that at this point ms_bscat_double contains the
	   integral in Eloranta's formulation */
	ms_bscat_air_double[i] = ms_bscat_air_single[i] * fov_factor 
	  * (fs_optical_depth - ms_bscat_double[i]*drange);
	/* Now set ms_bscat_double to its proper value */
	ms_bscat_double[i] = afactor
	  * ms_bscat_single[i] * fov_factor 
	  * (fs_optical_depth - ms_bscat_double[i]*drange);
      }
      else if (ms_E_once[i] > 0.0) {
	/* Inexact double scattering formulation */
	ms_real double_ratio = 2.0 * total_range2 * mu2_receiver * fov_factor 
	  * INT_TOPHAT(ms_E_once[i], 1, ms_Emu2_once[i],
		       mu2_receiver*total_range2);
	ms_bscat_double[i] = ms_bscat_single[i] * double_ratio;
	ms_bscat_air_double[i] = ms_bscat_air_single[i] * double_ratio;
      }

      /* Now work out the backscatter at gate i */
      if (ms_E_multi[i] > 0.0) {
	if (config->options & MS_APPROXIMATE_EXPONENTIAL) {
	  integral_ratio = 2.0 * total_range2 * mu2_receiver * fov_factor * (
		      INT_GAUSSIAN(ms_E_multi[i], 1,
				 ms_Emu2_multi[i], mu2_receiver*total_range2));
	}
	else {
	  integral_ratio = afactor
	    * 2.0 * total_range2 * mu2_receiver * fov_factor * (
		      INT_TOPHAT(ms_E_multi[i], 1,
				 ms_Emu2_multi[i], mu2_receiver*total_range2));
	}
      }
      if (config->options & MS_DOUBLE_SCATTERING_ONLY) {
	if (bscat_air_out) {
	  bscat_out[i] = ms_bscat_single[i] + ms_bscat_double[i];
	  bscat_air_out[i] = ms_bscat_air_single[i] + ms_bscat_air_double[i];
	}
	else {
	  bscat_out[i] = ms_bscat_single[i] + ms_bscat_air_single[i]
	    + ms_bscat_double[i] + ms_bscat_air_double[i];
	}
      }
      else {
	if (bscat_air_out) {
	  bscat_out[i] = ms_bscat_double[i] 
	    + ms_bscat_single[i] * (1.0 + integral_ratio);
	  bscat_air_out[i] = ms_bscat_air_double[i] 
	    + ms_bscat_air_single[i] * (1.0 + integral_ratio);
	}
	else {
	  bscat_out[i] = ms_bscat_double[i] + ms_bscat_air_double[i]
	    + (ms_bscat_single[i]+ms_bscat_air_single[i]) 
	     * (1.0 + integral_ratio);
	}
	/* Save the individual backscatter components */
	ms_bscat_multi[i]  = ms_bscat_single[i]*integral_ratio;
	ms_bscat_air_multi[i]  = ms_bscat_air_single[i]*integral_ratio;
      }
    }
    else {
      bscat_out[i] = 0.0;
      if (bscat_air_out) {
	bscat_air_out[i] = 0.0;
      }
      /* In case still set by Eloranta's double-scattering method */
      ms_bscat_double[i] = 0.0;
    }
  }
  return MS_SUCCESS;
}
