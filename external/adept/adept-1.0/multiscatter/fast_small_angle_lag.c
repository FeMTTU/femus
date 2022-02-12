/* fast_small_angle_lag.c -- Fast small-angle algorithm with pulse lag

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

/* Use the photon variance-covariance approach to calculate the
   backscatter with O(N) efficiency */
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
    ms_real *lag_out)         /* Pulse lag, m */
{
  int i;
  ms_real rho_transmitter2 = instrument.rho_transmitter
    *instrument.rho_transmitter;
  ms_real bscat_ext_ratio_air = 3.0 / (8.0 * MS_PI);
  ms_real rho_factor;

  /* Properties at previous half-levels */
  ms_real midRange_prev = fabs(instrument.altitude-1.5*range[0]+0.5*range[1]);
  ms_real M_prev = 0.0;
  ms_real afactor_prev = 1.0;
  /* Initial properties of the total distribution */
  ms_real width2 = midRange_prev*midRange_prev*rho_transmitter2;
  ms_real zeta2 = rho_transmitter2;
  ms_real cov = midRange_prev*rho_transmitter2;
  ms_real trans = 1.0;
  /* Initial properties of the "reduced" total distribution */
  ms_real width2r = width2;
  ms_real zeta2r = zeta2;
  ms_real covr = cov;
  ms_real transr = trans;
  /* Initial properties of the unscattered distribution */
  ms_real transu = trans;
  ms_real transu_prev = transu;

  if (instrument.receiver_type == MS_TOP_HAT) {
    rho_factor = 1.0/(1.0-exp(-instrument.rho_receiver[0]
			      *instrument.rho_receiver[0]/
			      rho_transmitter2));
  }
  else {
    rho_factor = 1.0 + rho_transmitter2
      / (instrument.rho_receiver[0]*instrument.rho_receiver[0]);
  }

  /* Main loop (note that this is the only one!) */
  for (i = 0; i < n; i++) {
    /* Variables at full levels */
    ms_real drange = ms_get_drange(i, n, range);
    ms_real drange2 = drange*drange;
    ms_real drange3 = drange2*drange;
    ms_real Theta = instrument.wavelength/(MS_PI*radius[i]);
    ms_real Theta2 = Theta*Theta;
    ms_real bscat_unattenuated = ext[i]/ext_bscat_ratio[i];
    ms_real bscat_air_unattenuated = ext_air[i]*bscat_ext_ratio_air;
    ms_real layer_od = 2.0*(ext[i]+ext_air[i])*drange;
    /* Variables at half levels */
    ms_real Range = fabs(instrument.altitude-range[i]);
    ms_real midRange = Range+drange*0.5;
    ms_real midRange2 = midRange*midRange;
    ms_real width2a_max = midRange2
      *instrument.rho_receiver[0]*instrument.rho_receiver[0];
    ms_real M = 0.0, transa;
    ms_real afactor = 1.0, afactor_next = 1.0;

    /* Properties of the unscattered distribution */
    ms_real width2u = midRange2*rho_transmitter2;
    ms_real covu = midRange*rho_transmitter2;
    
    /* Step the unscattered transmittance forward */
    transu = transu*exp(-layer_od);

    /* Step the total distribution at half-levels forward */
    width2 = width2 + zeta2*drange2
      + 2.0*cov*drange + ext[i]*Theta2*drange3*0.333333;
    // ERROR IN ORIGINAL CODE
    //    cov = cov + zeta2*drange + ext[i]*Theta2*drange*0.5;
    cov = cov + zeta2*drange + ext[i]*Theta2*drange2*0.5;
    zeta2 = zeta2 + ext[i]*drange*Theta2;
    trans = trans*exp(-(ext[i]+2.0*ext_air[i])*drange);
    
    /* Step the "reduced" total distribution at half-levels
       forward */
    width2r = width2r + zeta2r*drange2
      + 2.0*covr*drange + ext[i]*Theta2*drange3*0.333333;
    // ERROR IN ORIGINAL CODE
    //    covr = covr + zeta2r*drange + ext[i]*Theta2*drange*0.5;
    covr = covr + zeta2r*drange + ext[i]*Theta2*drange2*0.5;
    zeta2r = zeta2r + ext[i]*drange*Theta2;
    transr = transr*exp(-(ext[i]+2.0*ext_air[i])*drange);

    /* Calculate the properties of the first Gaussian */
    transa = transr-transu;
    if (transa < 0.0) {
      transa = 0.0;
    }
    if (transa > transr*1.0e-12) {
      /* There has been some forward scattering */
      int is_anisotropic = droplet_fraction && pristine_ice_fraction 
	&& (droplet_fraction[i] > 0.0 || pristine_ice_fraction[i] > 0.0);
      int is_anisotropic_next = i < n-1 && droplet_fraction && pristine_ice_fraction
	&& (droplet_fraction[i+1] > 0.0 || pristine_ice_fraction[i+1] > 0.0);
      ms_real width2a = (transr*width2r - transu*width2u)/transa;
      ms_real zeta2a = (transr*zeta2r - transu*rho_transmitter2)/transa;
      ms_real cova = (transr*covr - transu*covu)/transa;
      ms_real transb = trans - transr;
      ms_real bscat_factor;
      if (width2a > width2a_max) {
	/* The first Gaussian is larger than the receiver FOV and
	   needs adjusting */
	ms_real correlation2a = cova*cova/(width2a*zeta2a);
	ms_real factor = width2a_max/width2a;
	width2a = width2a_max;
	transa = transa*factor;
	cova = cova*factor;
	zeta2a = zeta2a*(correlation2a*factor + 1-correlation2a);
	/* Recalculate "reduced" total distribution */
	transr = transu+transa;
	width2r = (transa*width2a + transu*width2u)/transr;
	zeta2r = (transa*zeta2a + transu*rho_transmitter2)/transr;
	covr = (transa*cova + transu*covu)/transr;
      }
      if (is_anisotropic || is_anisotropic_next) {
	/* Calculate anisotropic scaling factor */
	/* First calculate the properties of all the forward scattered
	   photons together */
	ms_real transg = trans-transu;
	ms_real width2g = (trans*width2 - transu*width2u)/transg;
	ms_real zeta2g = (trans*zeta2 - transu*rho_transmitter2)/transg;
	ms_real covg = (trans*cov - transu*covu)/transg;
	if (is_anisotropic) {
	  afactor = ms_anisotropic_factor(width2g, zeta2g, covg,
				       Theta2, midRange, width2a_max,
				       bscat_unattenuated,
				       bscat_air_unattenuated,
				       droplet_fraction[i], pristine_ice_fraction[i]);
	}
	if (is_anisotropic_next) {
	  /* The backscatter for the next full level needs an
	     anisotropic factor calculated at the current half level
	     so it is easiest to do now when the variables describing
	     the distribution are available */
	  ms_real Theta_next = instrument.wavelength/(MS_PI*radius[i+1]);
	  afactor_next = ms_anisotropic_factor(width2g, zeta2g, covg,
					    Theta_next*Theta_next,
					    midRange, width2a_max,
					    ext[i+1]/ext_bscat_ratio[i],
					    ext_air[i+1]*bscat_ext_ratio_air,
					    droplet_fraction[i+1],
					    pristine_ice_fraction[i+1]);
	}
      }

      if (transb > 0.0) {
	/* We have two Gaussians */
	ms_real width2b = (trans*width2 - transr*width2r)/transb;

	/* Calculate M, the isotropic-equivalent backscatter
	   enhancement factor at the current half-gate */
	if (instrument.receiver_type == MS_TOP_HAT) {
	  M = ((1.0-exp(-width2a_max/width2a))*transa
	       +(1.0-exp(-width2a_max/width2b))*transb)
	    * rho_factor / transu;
	}
	else {
	  M = rho_factor * (transa / (1.0 + width2a/width2a_max)
			    +transb / (1.0+width2b/width2a_max)) / transu;
	}
      }
      /* Otherwise we have one Gaussian */
      else if (instrument.receiver_type == MS_TOP_HAT) {
	M = rho_factor*(1.0-exp(-width2a_max/width2a))*transa/transu;
      }
      else {
	M = rho_factor*transa/((1.0 + width2a/width2a_max)*transu);
      }

      /* The following operation calculates the apparent backscatter
	 given the "backscatter enhancement factor" at the previous
	 half-gate (afactor_prev*M_prev) and the same at the current
	 half-gate (afactor*M). This is done by assuming that the
	 variable (1+afactor*M) varies exponentially within the gate,
	 while the unattenuated backscatter coefficient remains
	 constant. */
      bscat_factor = transu_prev*(1.0+afactor_prev*M_prev
				  - (1+afactor*M)*exp(-layer_od))
	/(layer_od-log((1+afactor*M)/(1+afactor_prev*M_prev)));
      
      if (bscat_air_out) {
	bscat_out[i] = bscat_factor*bscat_unattenuated;
	bscat_air_out[i] = bscat_factor*bscat_air_unattenuated;
      }
      else {
	bscat_out[i] = bscat_factor
	  * (bscat_unattenuated+bscat_air_unattenuated);
      }

      if (lag_out) {
	if (transb <= 0.0) {
	  lag_out[i] = cov - covu;
	}
	else {
	  ms_real width2b = (trans*width2 - transr*width2r)/transb;
	  ms_real covb = (trans*cov - transr*covr)/transb;
	  lag_out[i] = (covu*transu + cova*transa
			+ covb*transb*width2a_max/width2b
			*width2a_max/width2b)
	    / (transu + transa + transb*width2a_max/width2b) - covu;
	}
      }
    }
    else if (layer_od > 0.0) {
      /* No forward scattering but some attenuation */
      ms_real bscat_factor = transu_prev*(1.0-exp(-layer_od))/layer_od;
      if (bscat_air_out) {
	bscat_out[i] = bscat_factor*bscat_unattenuated;
	bscat_air_out[i] = bscat_factor*bscat_air_unattenuated;
      }
      else {
	bscat_out[i] = bscat_factor
	  * (bscat_unattenuated+bscat_air_unattenuated);
      }
      if (lag_out) {
	lag_out[i] = 0.0;
      }
    }
    else {
      bscat_out[i] = 0.0;
      if (bscat_air_out) {
	bscat_air_out[i] = 0.0;
      }
      if (lag_out) {
	lag_out[i] = 0.0;
      }
    }
    /* Set variables for the previous half-level */
    midRange_prev = midRange;
    transu_prev = transu;
    M_prev = M;
    afactor_prev = afactor_next;
  }

  return MS_SUCCESS;
}

