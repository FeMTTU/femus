/* singlescatter_AD.c -- Adjoint of single scattering calculation

   Copyright (C) 2009-2011 Robin Hogan <r.j.hogan@reading.ac.uk> 

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

#include <math.h>

#include "ms.h"

/* Perform single scattering calculation with no delta scaling. */
int
ms_singlescatter_AD(
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
    ms_real *ext_bscat_ratio_AD)
{
  ms_real transmittance[n+1];
  ms_real optical_depth = 0.0;
  ms_real bscat_ext_ratio_air = 3.0 / (8.0 * MS_PI);
  int i;

  /* Assume that for wavelengths longer than 1 micron the gaseous
     extinction is entirely due to absorption rather than Rayleigh
     scattering */
  if (instrument.wavelength > 1e-6) {
    bscat_ext_ratio_air = 0.0;
  }

  transmittance[0] = 1.0;

  /* *** FORWARD CALCULATION *** */
  if (ext_air) {
    /* Scattering atmospheric gases are present */
    for (i = 0; i < n; i++) {
      /* optical depth of layer */
      ms_real two_drange = 2.0*ms_get_drange(i, n, range);
      ms_real optical_depth_layer = (ext_air[i]+ext[i])*two_drange;
      ms_real factor;
      if (optical_depth_layer > 0.0) {
	optical_depth += optical_depth_layer;
	transmittance[i+1] = exp(-optical_depth);
	factor = (transmittance[i]-transmittance[i+1])/optical_depth_layer;
      }
      else {
	factor = transmittance[i+1] = transmittance[i];
      }

      if (bscat_air_out) {
	/* Separate particulate and molecular returns */
	bscat_out[i] = factor*ext[i]/ext_bscat_ratio[i];
	bscat_air_out[i] = factor*ext_air[i]*bscat_ext_ratio_air;
      }
      else {
	bscat_out[i] = factor*(ext[i]/ext_bscat_ratio[i] 
			       + ext_air[i]*bscat_ext_ratio_air);
      }
    }
  }
  else {
    /* No atmospheric gases to consider for scattering, although they
       might attenuate */
    for (i = 0; i < n; i++) {
      ms_real two_drange = 2.0*ms_get_drange(i, n, range);
      ms_real optical_depth_layer;
      if (ext_air) {
	optical_depth_layer = (ext_air[i]+ext[i])*two_drange;
      }
      else {
	optical_depth_layer = ext[i]*two_drange;
      }

      if (optical_depth_layer > 0.0) {
	optical_depth += optical_depth_layer;
	transmittance[i+1] = exp(-optical_depth);
	bscat_out[i] = (transmittance[i]-transmittance[i+1])*ext[i]
	  / (optical_depth_layer*ext_bscat_ratio[i]);
      }
      else {
	bscat_out[i] = 0.0;
	transmittance[i+1] = transmittance[i];
      }
    }
    if (bscat_air_out) {
      for (i = 0; i < n; i++) {
	bscat_air_out[i] = 0.0;
      }
    }
  }

    /* *** ADJOINT CALCULATION *** */
  if (ext_air) {
    /* Scattering atmospheric gases are present */
    ms_real sum = 0.0;
    for (i = n-1; i >= 0; i--) {
      ms_real two_drange = 2.0*ms_get_drange(i, n, range);
      if (bscat_air_out) {
	/* Particulate and air backscatters are being separated */
	ms_real B, Bair = 0.0;
	if (ext[i] > 0.0) {
	  B = bscat_out[i]/ext[i]
	    + (ext[i]*transmittance[i+1]/ext_bscat_ratio[i] - bscat_out[i])
	    / (ext[i] + ext_air[i]);
	}
	else if (i > 0) {
	  B = transmittance[i]/ext_bscat_ratio[i];
	}
	else {
	  B = 1.0/ext_bscat_ratio[i];
	}


	if (ext_air[i] > 0.0) {
	  Bair = (ext_air[i]*bscat_ext_ratio_air*transmittance[i+1]
		  -bscat_air_out[i]) / (ext[i] + ext_air[i]);
	}

	ext_AD[i] += B*bscat_AD[i] + Bair*bscat_air_AD[i] - sum;
	sum += two_drange*(bscat_out[i]*bscat_AD[i]
			   + bscat_air_out[i]*bscat_air_AD[i]);
	ext_bscat_ratio_AD[i] -= bscat_AD[i]*bscat_out[i]/ext_bscat_ratio[i];
      }
      else {
	/* "bscat_out" and "bscat_AD" refer to the total backscatter */
	ms_real B;
	if (ext[i] > 0.0) {
	  ms_real bscat_particulate_fraction = 1.0
	    /(1.0+bscat_ext_ratio_air*ext_air[i]
	      *ext_bscat_ratio[i]/ext[i]);
	  ext_bscat_ratio_AD[i] -= bscat_AD[i]*bscat_out[i]
	    *bscat_particulate_fraction/ext_bscat_ratio[i];
	  B = bscat_particulate_fraction*bscat_out[i]/ext[i]
	    + (transmittance[i+1]
	       *(ext[i]/ext_bscat_ratio[i]
		 +bscat_ext_ratio_air*ext_air[i]) - bscat_out[i])
	    / (ext[i] + ext_air[i]);
	}
	else if (i > 0) {
	  B = transmittance[i-1]/ext_bscat_ratio[i];
	}
	else {
	  B = 1.0/ext_bscat_ratio[i];
	}

	ext_AD[i] += B*bscat_AD[i] - sum;
	sum += two_drange*(bscat_out[i]*bscat_AD[i]);
      }
    }
  }
  else {
    /* No atmospheric gases to consider for scattering, although they
       might attenuate */ 
    ms_real two_drange = 2.0*ms_get_drange(i, n, range);
    ms_real sum = 0.0;
    for (i = n-1; i >= 0; i--) {
      if (ext[i] > 0.0) {
	ms_real B;
	if (ext_air) {
	  B = bscat_out[i]/ext[i]
	    + (ext[i]*transmittance[i+1]/ext_bscat_ratio[i] - bscat_out[i])
	    / (ext[i] + ext_air[i]);
	}
	else {
	  B = transmittance[i+1]/ext_bscat_ratio[i];
	}
	ext_AD[i] += B*bscat_AD[i] - sum;
	sum += two_drange*bscat_out[i]*bscat_AD[i];
	ext_bscat_ratio_AD[i] -= bscat_AD[i]*bscat_out[i]/ext_bscat_ratio[i];
      }
    }
  }

  return MS_SUCCESS;
}


