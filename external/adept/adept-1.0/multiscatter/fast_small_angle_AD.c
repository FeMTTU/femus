/* fast_small_angle_AD.c -- Adjoint of fast small-angle multiple scattering algorithm

   Copyright (C) 2007-2011 Robin Hogan <r.j.hogan@reading.ac.uk> 

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


/* This is the adjoint of the photon variance-covariance algorithm to
   calculate the backscatter with O(N) efficiency. To understand what
   the original code does, read:

 Hogan, R. J., 2008: Fast lidar and radar multiple-scattering models -
 1. Small-angle scattering using the photon variance-covariance
 method.  J. Atmos. Sci., 65, 3621-3635.

*/

int
ms_fast_small_angle_AD(
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
    ms_real *pristine_ice_fraction_AD)
{
  int i, n1 = n+1;
  ms_real rho_transmitter2 = instrument.rho_transmitter
    *instrument.rho_transmitter;
  ms_real bscat_ext_ratio_air = 3.0 / (8.0 * MS_PI);
  ms_real rho_factor;

  /* Properties at previous half-levels */
  ms_real midRange_prev = fabs(instrument.altitude-1.5*range[0]+0.5*range[1]);

  ms_real width2[n1];
  ms_real zeta2[n1];
  ms_real cov[n1];
  ms_real trans[n1];
  ms_real width2r[n1];
  ms_real zeta2r[n1];
  ms_real covr[n1];
  ms_real transr[n1];
  ms_real width2a[n1];
  ms_real zeta2a[n1];
  ms_real cova[n1];
  ms_real transa[n1];
  ms_real transu[n1];
  ms_real transb[n1];
  ms_real M[n1];
  ms_real afactor_far[n]; /* "afactor" in the non-adjoint version */
  ms_real afactor_near[n];/* afactor_prev and afactor_next in non-adj */

  ms_real width2r_orig[n1];
  ms_real zeta2r_orig[n1];
  ms_real covr_orig[n1];
  ms_real transr_orig[n1];
  ms_real width2a_orig[n1];
  ms_real zeta2a_orig[n1];
  ms_real cova_orig[n1];
  ms_real transa_orig[n1];

  ms_real width2_AD[n1];
  ms_real zeta2_AD[n1];
  ms_real cov_AD[n1];
  ms_real trans_AD[n1];
  ms_real width2r_AD[n1];
  ms_real zeta2r_AD[n1];
  ms_real covr_AD[n1];
  ms_real transr_AD[n1];
  ms_real width2a_AD[n1];
  ms_real zeta2a_AD[n1];
  ms_real cova_AD[n1];
  ms_real transa_AD[n1];
  ms_real transu_AD[n1];
  ms_real transb_AD[n1];
  ms_real M_AD[n1];
  ms_real afactor_far_AD[n];
  ms_real afactor_near_AD[n];

  ms_real width2r_orig_AD[n1];
  ms_real zeta2r_orig_AD[n1];
  ms_real covr_orig_AD[n1];
  ms_real transr_orig_AD[n1];
  ms_real width2a_orig_AD[n1];
  ms_real zeta2a_orig_AD[n1];
  ms_real cova_orig_AD[n1];
  ms_real transa_orig_AD[n1];

  /* Assume that for wavelengths longer than 1 micron the gaseous
     extinction is entirely due to absorption rather than Rayleigh
     scattering */
  if (instrument.wavelength > 1e-6) {
    bscat_ext_ratio_air = 0.0;
  }

  for (i = 0; i < n1; i++) {
    width2_AD[i] = zeta2_AD[i] = cov_AD[i] = trans_AD[i] = 0.0;
    width2r_AD[i] = zeta2r_AD[i] = covr_AD[i] = transr_AD[i] = 0.0;
    width2a_AD[i] = zeta2a_AD[i] = cova_AD[i] = transa_AD[i] = 0.0;
    width2r_orig_AD[i] = zeta2r_orig_AD[i]
      = covr_orig_AD[i] = transr_orig_AD[i] = 0.0;
    width2a_orig_AD[i] = zeta2a_orig_AD[i]
      = cova_orig_AD[i] = transa_orig_AD[i] = 0.0;
    transu_AD[i] = transb_AD[i] = M_AD[i] = 0.0;
  }
  for (i = 0; i < n; i++) {
    afactor_far_AD[i] = afactor_near_AD[i] = 0.0;
  }

  /* Initial properties of the total distribution */
  width2[0] = midRange_prev*midRange_prev*rho_transmitter2;
  zeta2[0] = rho_transmitter2;
  cov[0] = midRange_prev*rho_transmitter2;
  trans[0] = 1.0;
  /* Initial properties of the "reduced" total distribution */
  width2r[0] = width2r_orig[0] = width2[0];
  zeta2r[0] = zeta2r_orig[0] = zeta2[0];
  covr[0] = covr_orig[0] = cov[0];
  transr[0] = transr_orig[0] = trans[0];
  /* Initial properties of the unscattered distribution */
  transu[0] = trans[0];

  M[0] = 0.0;
  afactor_far[0] = 1.0;
  afactor_near[0] = 1.0;

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
    // SECTION 0
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

    /* Properties of the unscattered distribution */
    ms_real width2u = midRange2*rho_transmitter2;
    ms_real covu = midRange*rho_transmitter2;

    int i1 = i+1;

    // SECTION 1
    M[i1] = 0.0;
    afactor_far[i] = 1.0;
    afactor_near[i1] = 1.0;
    
    /* Step the unscattered transmittance forward */
    transu[i1] = transu[i]*exp(-layer_od);

    /* Step the total distribution at half-levels forward */
    width2[i1] = width2[i] + zeta2[i]*drange2
      + 2.0*cov[i]*drange + ext[i]*Theta2*drange3*0.333333;
    // ERROR IN ORIGINAL CODE!
    //    cov[i1] = cov[i] + zeta2[i]*drange + ext[i]*Theta2*drange*0.5;
    cov[i1] = cov[i] + zeta2[i]*drange + ext[i]*Theta2*drange2*0.5;
    zeta2[i1] = zeta2[i] + ext[i]*drange*Theta2;
    trans[i1] = trans[i]*exp(-(ext[i]+2.0*ext_air[i])*drange);
    
    /* Step the "reduced" total distribution at half-levels
       forward */
    width2r_orig[i1] = width2r[i] + zeta2r[i]*drange2
      + 2.0*covr[i]*drange + ext[i]*Theta2*drange3*0.333333;
    // ERROR IN ORIGINAL CODE!
    //    covr_orig[i1] = covr[i] + zeta2r[i]*drange + ext[i]*Theta2*drange*0.5;
    covr_orig[i1] = covr[i] + zeta2r[i]*drange + ext[i]*Theta2*drange2*0.5;
    zeta2r_orig[i1] = zeta2r[i] + ext[i]*drange*Theta2;
    transr_orig[i1] = transr[i]*exp(-(ext[i]+2.0*ext_air[i])*drange);

    /* Calculate the properties of the first Gaussian */
    transa_orig[i1] = transr_orig[i1]-transu[i1];
    if (transa_orig[i1] < 0.0) {
      transa_orig[i1] = 0.0;
    }

    if (transa_orig[i1] > transr_orig[i1]*1.0e-12) {
      // SECTION 2
      // SECTION 2.1
      /* There has been some forward scattering */
      int is_anisotropic = droplet_fraction && pristine_ice_fraction
	&& (droplet_fraction[i] > 0.0 || pristine_ice_fraction[i] > 0.0);
      int is_anisotropic_next
	= i < n-1 && droplet_fraction && pristine_ice_fraction
	&& (droplet_fraction[i1] > 0.0 || pristine_ice_fraction[i1] > 0.0);
      width2a_orig[i1] = (transr_orig[i1]*width2r_orig[i1]
			  - transu[i1]*width2u)/transa_orig[i1];
      zeta2a_orig[i1] = (transr_orig[i1]*zeta2r_orig[i1]
			 - transu[i1]*rho_transmitter2)/transa_orig[i1];
      cova_orig[i1] = (transr_orig[i1]*covr_orig[i1]
		       - transu[i1]*covu)/transa_orig[i1];
      transb[i1] = trans[i1] - transr_orig[i1];
      // SECTION 2.2
      if (width2a_orig[i1] > width2a_max) {
	// SECTION 2.2.1
	/* The first Gaussian is larger than the receiver FOV and
	   needs adjusting */
	ms_real correlation2a = cova_orig[i1]*cova_orig[i1]
	  /(width2a_orig[i1]*zeta2a_orig[i1]);
	ms_real factor = width2a_max/width2a_orig[i1];
	width2a[i1] = width2a_max;
	transa[i1] = transa_orig[i1]*factor;
	cova[i1] = cova_orig[i1]*factor;
	zeta2a[i1] = zeta2a_orig[i1]*(correlation2a*factor
				      + 1.0-correlation2a);
	/* Recalculate "reduced" total distribution */
	transr[i1] = transu[i1]+transa[i1];
	width2r[i1] = (transa[i1]*width2a[i1]
			+ transu[i1]*width2u)/transr[i1];
	zeta2r[i1] = (transa[i1]*zeta2a[i1]
		       + transu[i1]*rho_transmitter2)/transr[i1];
	covr[i1] = (transa[i1]*cova[i1] + transu[i1]*covu)/transr[i1];
      }
      else {
	// SECTION 2.2.2
	transa[i1] = transa_orig[i1];
	width2a[i1] = width2a_orig[i1];
	cova[i1] = cova_orig[i1];
	zeta2a[i1] = zeta2a_orig[i1];

	transr[i1] = transr_orig[i1];
	width2r[i1] = width2r_orig[i1];
	covr[i1] = covr_orig[i1];
	zeta2r[i1] = zeta2r_orig[i1];
      }
      // SECTION 2.3
      if (is_anisotropic || is_anisotropic_next) {
	// SECTION 2.3.1
	/* Calculate anisotropic scaling factor */
	/* First calculate the properties of all the forward scattered
	   photons together */
	ms_real transg = trans[i1]-transu[i1];
	ms_real width2g = (trans[i1]*width2[i1]
			   - transu[i1]*width2u)/transg;
	ms_real zeta2g = (trans[i1]*zeta2[i1]
			  - transu[i1]*rho_transmitter2)/transg;
	ms_real covg = (trans[i1]*cov[i1] - transu[i1]*covu)/transg;
	// SECTION 2.3.2
	if (is_anisotropic) {
	  afactor_far[i] = ms_anisotropic_factor(width2g, zeta2g, covg,
					  Theta2, midRange, width2a_max,
					  bscat_unattenuated,
					  bscat_air_unattenuated,
					  droplet_fraction[i],
					  pristine_ice_fraction[i]);
	}
	if (is_anisotropic_next) {
	  /* The backscatter for the next full level needs an
	     anisotropic factor calculated at the current half level
	     so it is easiest to do now when the variables describing
	     the distribution are available */
	  ms_real Theta_next = instrument.wavelength/(MS_PI*radius[i1]);
	  afactor_near[i1] = ms_anisotropic_factor(width2g, zeta2g, covg,
					       Theta_next*Theta_next,
					       midRange, width2a_max,
					       ext[i1]/ext_bscat_ratio[i1],
					       ext_air[i1]*bscat_ext_ratio_air,
					       droplet_fraction[i1],
					       pristine_ice_fraction[i1]);
	}
      }

      // SECTION 2.4
      if (transb[i1] > 0.0) {
	// SECTION 2.4.1
	/* We have two Gaussians */
	ms_real width2b = (trans[i1]*width2[i1]
			   - transr[i1]*width2r[i1])/transb[i1];
	/* Sometimes width2b can be 0 if transb is very small */
	if (width2b <= width2a[i1]) {
	  width2b = width2a[i1];
	}
	/* Calculate M, the isotropic-equivalent backscatter
	   enhancement factor at the current half-gate */
	if (instrument.receiver_type == MS_TOP_HAT) {
	  M[i1] = ((1.0-exp(-width2a_max/width2a[i1]))*transa[i1]
	       +(1.0-exp(-width2a_max/width2b))*transb[i1])
	    * rho_factor / transu[i1];
	}
	else {
	  M[i1] = rho_factor * (transa[i1] / (1.0 + width2a[i1]/width2a_max)
				+transb[i1] / (1.0 + width2b/width2a_max))
	    / transu[i1];
	}
      }
      /* Otherwise we have one Gaussian */
      else if (instrument.receiver_type == MS_TOP_HAT) {
	// SECTION 2.4.2
	M[i1] = rho_factor*(1.0-exp(-width2a_max/width2a[i1]))
	  *transa[i1]/transu[i1];
      }
      else {
	// SECTION 2.4.3
	M[i1] = rho_factor*transa[i1]
	  /((1.0 + width2a[i1]/width2a_max)*transu[i1]);
      }
    
      // SECTION 2.5
      /* The following operation calculates the apparent backscatter
	 given the "backscatter enhancement factor" at the previous
	 half-gate (afactor_near[i]*M[i]) and the same at the current
	 half-gate (afactor_far[i]*M[i1]). This is done by assuming
	 that the variable (1+afactor*M) varies exponentially within
	 the gate, while the unattenuated backscatter coefficient
	 remains constant. */
      ms_real bscat_factor = transu[i]*(1.0+afactor_near[i]*M[i]
				- (1.0+afactor_far[i]*M[i1])*exp(-layer_od))
	/(layer_od-log((1.0+afactor_far[i]*M[i1])/(1.0+afactor_near[i]*M[i])));
      
      if (bscat_air_out) {
	bscat_out[i] = bscat_factor*bscat_unattenuated;
	bscat_air_out[i] = bscat_factor*bscat_air_unattenuated;
      }
      else {
	bscat_out[i] = bscat_factor
	  * (bscat_unattenuated+bscat_air_unattenuated);
      }
    }
    else {
      if (layer_od > 0.0) {
	// SECTION 3
	/* No forward scattering but some attenuation */
	ms_real bscat_factor = transu[i]*(1.0-exp(-layer_od))/layer_od;
	if (bscat_air_out) {
	  bscat_out[i] = bscat_factor*bscat_unattenuated;
	  bscat_air_out[i] = bscat_factor*bscat_air_unattenuated;
	}
	else {
	  bscat_out[i] = bscat_factor
	    * (bscat_unattenuated+bscat_air_unattenuated);
	}
      }
      else {
	// SECTION 4
	bscat_out[i] = 0.0;
	if (bscat_air_out) {
	  bscat_air_out[i] = 0.0;
	}
      }
      // SECTION 5
      transr[i1] = transr_orig[i1];
      width2r[i1] = width2r_orig[i1];
      covr[i1] = covr_orig[i1];
      zeta2r[i1] = zeta2r_orig[i1];
    }
  }

  /* *** ADJOINT CALCULATION *** */

  for (i = n-1; i >= 0; i--) {
    /* Variables at full levels */
    ms_real drange = ms_get_drange(i, n, range);
    ms_real drange2 = drange*drange;
    ms_real drange3 = drange2*drange;
    ms_real Theta = instrument.wavelength/(MS_PI*radius[i]);
    ms_real Theta2 = Theta*Theta;
    ms_real bscat_unattenuated = ext[i]/ext_bscat_ratio[i];
    ms_real bscat_air_unattenuated = ext_air[i]*bscat_ext_ratio_air;
    ms_real layer_od = 2.0*(ext[i]+ext_air[i])*drange;
    ms_real inv_exp_layer_od = exp(-layer_od);
    /* Variables at half levels */
    ms_real Range = fabs(instrument.altitude-range[i]);
    ms_real midRange = Range+drange*0.5;
    ms_real midRange2 = midRange*midRange;
    ms_real width2a_max = midRange2
      *instrument.rho_receiver[0]*instrument.rho_receiver[0];

    /* Properties of the unscattered distribution */
    ms_real width2u = midRange2*rho_transmitter2;
    ms_real covu = midRange*rho_transmitter2;

    ms_real Theta2_AD = 0.0;
    ms_real layer_od_AD = 0.0;
    ms_real bscat_unattenuated_AD = 0.0;
    ms_real bscat_air_unattenuated_AD = 0.0;


    int i1 = i+1;

    if (transa_orig[i1] > transr[i1]*1.0e-12) {
      // SECTION 2
      /* There has been some forward scattering */
      int is_anisotropic = droplet_fraction && pristine_ice_fraction 
	&& (droplet_fraction[i] > 0.0 || pristine_ice_fraction[i] > 0.0);
      int is_anisotropic_next 
	= i < n-1 && droplet_fraction && pristine_ice_fraction
	&& (droplet_fraction[i1] > 0.0 || pristine_ice_fraction[i1] > 0.0);
      ms_real num = transu[i]*(1.0+afactor_near[i]*M[i]
			       - (1.0+afactor_far[i]*M[i1])
			       *inv_exp_layer_od);
      ms_real denom = layer_od-log((1.0+afactor_far[i]*M[i1])
				   /(1.0+afactor_near[i]*M[i]));
      ms_real bscat_factor = num/denom;
      ms_real bscat_factor_AD = 0.0;

      // SECTION 2.5
      /* The following operation calculates the apparent backscatter
	 given the "backscatter enhancement factor" at the previous
	 half-gate (afactor_near[i]*M[i]) and the same at the current
	 half-gate (afactor_far[i]*M[i1]). This is done by assuming that the
	 variable (1+afactor*M) varies exponentially within the gate,
	 while the unattenuated backscatter coefficient remains
	 constant. */
      
      if (bscat_air_out) {
	bscat_unattenuated_AD += bscat_factor*bscat_AD[i];
	bscat_factor_AD += bscat_unattenuated*bscat_AD[i];
	/* Note that because bscat_AD will be used by the wide-angle
	   part of the code, we cannot set it to zero after use;
	   likewise with bscat_air_AD */
	//	bscat_AD[i] = 0.0;

	bscat_air_unattenuated_AD += bscat_factor*bscat_air_AD[i];
	bscat_factor_AD += bscat_air_unattenuated*bscat_air_AD[i];
	//	bscat_air_AD[i] = 0.0;
      }
      else {
	bscat_unattenuated_AD += bscat_factor*bscat_AD[i];
	bscat_air_unattenuated_AD += bscat_factor*bscat_AD[i];
	bscat_factor_AD += (bscat_unattenuated
			    +bscat_air_unattenuated)*bscat_AD[i];
	//	bscat_AD[i] = 0.0;
      }
    
      ms_real factor1 = (transu[i] - bscat_factor
			 /(1.0+afactor_near[i]*M[i]))/denom;
      ms_real factor2 = (-transu[i]*inv_exp_layer_od + bscat_factor
			 /(1.0+afactor_far[i]*M[i1]))/denom;
            
      afactor_near_AD[i] += bscat_factor_AD*M[i]*factor1;
      M_AD[i] += bscat_factor_AD*afactor_near[i]*factor1;
      afactor_far_AD[i] += bscat_factor_AD*M[i1]*factor2;
      M_AD[i1] += bscat_factor_AD*afactor_far[i]*factor2;
      layer_od_AD += bscat_factor_AD
	*(transu[i]*(1.0+afactor_far[i]*M[i1])*inv_exp_layer_od
	  -bscat_factor)/denom;
      transu_AD[i] += bscat_factor_AD*bscat_factor/transu[i];
      bscat_factor_AD = 0.0;
      
      // SECTION 2.4
      if (transb[i1] > 0.0) {
	// SECTION 2.4.1
	ms_real width2b = (trans[i1]*width2[i1]
			   - transr[i1]*width2r[i1])/transb[i1];
	ms_real width2b_AD = 0.0;
	if (width2b <= width2a[i1]) {
	  width2b = width2a[i1];
	}
	/* Calculate M, the isotropic-equivalent backscatter
	   enhancement factor at the current half-gate */
	if (instrument.receiver_type == MS_TOP_HAT) {
	  ms_real tmp_factor = rho_factor/transu[i1];
	  transa_AD[i1] += M_AD[i1]*tmp_factor
	    *(1.0-exp(-width2a_max/width2a[i1]));
	  width2a_AD[i1] -= M_AD[i1]*tmp_factor*transa[i1]
	    *exp(-width2a_max/width2a[i1])*width2a_max
	    /(width2a[i1]*width2a[i1]);
	  transb_AD[i1] += M_AD[i1]*tmp_factor
	    *(1.0-exp(-width2a_max/width2b));
	  width2b_AD -= M_AD[i1]*tmp_factor*transb[i1]
	    *exp(-width2a_max/width2b)*width2a_max
	    /(width2b*width2b);
	  transu_AD[i1] -= M_AD[i1]*M[i1]/transu[i1];
	  M_AD[i1] = 0.0;
	}
	else {
	  ms_real denom1 = 1.0+width2a[i1]/width2a_max;
	  ms_real denom2 = 1.0+width2b/width2a_max;
	  ms_real tmp_factor = rho_factor/transu[i1];
	  transa_AD[i1] += M_AD[i1]*tmp_factor/denom1;
	  transb_AD[i1] += M_AD[i1]*tmp_factor/denom2;
	  width2a_AD[i1] -= M_AD[i1]*tmp_factor*transa[i1]
	    /(denom1*denom1*width2a_max);
	  width2b_AD -= M_AD[i1]*tmp_factor*transb[i1]
	    /(denom2*denom2*width2a_max);
	  transu_AD[i1] -= M_AD[i1]*M[i1]/transu[i1];
	  M_AD[i1] = 0.0;
	}
	trans_AD[i1] += width2b_AD*width2[i1]/transb[i1];
	width2_AD[i1] += width2b_AD*trans[i1]/transb[i1];
	transr_AD[i1] -= width2b_AD*width2r[i1]/transb[i1];
	width2r_AD[i1] -= width2b_AD*transr[i1]/transb[i1];
	transb_AD[i1] -= width2b_AD*width2b/transb[i1];
	width2b_AD = 0.0;
      }
      /* Otherwise we have one Gaussian */
      else if (instrument.receiver_type == MS_TOP_HAT) {
	// SECTION 2.4.2
	transa_AD[i1] += M_AD[i1]*M[i1]/transa[i1];
	width2a_AD[i1] -= M_AD[i1]*rho_factor*transa[i1]
	  *exp(-width2a_max/width2a[i1])*width2a_max
	  /(width2a[i1]*width2a[i1]*transu[i1]);
	transu_AD[i1] -= M_AD[i1]*M[i1]/transu[i1];
	M_AD[i1] = 0.0;
      }
      else {
	// SECTION 2.4.3
	transa_AD[i1] += M_AD[i1]*M[i1]/transa[i1];
	width2a_AD[i1] -= M_AD[i1]*M[i1]
	  /((1.0+width2a[i1]/width2a_max)*width2a_max);
	transu_AD[i1] -= M_AD[i1]*M[i1]/transu[i1];
	M_AD[i1] = 0.0;
      }

      // SECTION 2.3
      if (is_anisotropic || is_anisotropic_next) {
	/* Calculate anisotropic scaling factor */
	/* First calculate the properties of all the forward scattered
	   photons together */
	ms_real transg = trans[i1]-transu[i1];
	ms_real width2g = (trans[i1]*width2[i1] 
			   - transu[i1]*width2u)/transg;
	ms_real zeta2g = (trans[i1]*zeta2[i1]
			  - transu[i1]*rho_transmitter2)/transg;
	ms_real covg = (trans[i1]*cov[i1] - transu[i1]*covu)/transg;
	ms_real width2g_AD = 0.0;
	ms_real zeta2g_AD = 0.0;
	ms_real covg_AD = 0.0;
	ms_real transg_AD = 0.0;
	// SECTION 2.3.2
	if (is_anisotropic) {
	  /* Ignore the return value from this function since we
	     already have calculated afactor_far[i] */
	  ms_anisotropic_factor_AD(width2g, zeta2g, covg,
				   Theta2, midRange, width2a_max,
				   bscat_unattenuated,
				   bscat_air_unattenuated,
				   droplet_fraction[i],
				   pristine_ice_fraction[i],
				   afactor_far_AD[i],
				   &width2g_AD,
				   &zeta2g_AD,
				   &covg_AD,
				   &Theta2_AD,
				   &bscat_unattenuated_AD,
				   &bscat_air_unattenuated_AD,
				   &(droplet_fraction_AD[i]),
				   &(pristine_ice_fraction_AD[i]));
	}
	if (is_anisotropic_next) {
	  // THIS SECTION NEEDS CHECKING
	  /* The backscatter for the next full level needs an
	     anisotropic factor calculated at the current half level
	     so it is easiest to do now when the variables describing
	     the distribution are available */
	  ms_real Theta_next = instrument.wavelength/(MS_PI*radius[i1]);
	  ms_real Theta2_next = Theta_next*Theta_next;
	  ms_real Theta2_next_AD = 0.0;
	  ms_real bscat_next_unattenuated = ext[i1]/ext_bscat_ratio[i1];
	  ms_real bscat_air_next_unattenuated
	    = ext_air[i1]*bscat_ext_ratio_air;
	  ms_real bscat_next_unattenuated_AD = 0.0;
	  ms_real bscat_air_next_unattenuated_AD = 0.0;
	  ms_anisotropic_factor_AD(width2g, zeta2g, covg,
				   Theta2_next,
				   midRange, width2a_max,
				   bscat_next_unattenuated,
				   bscat_air_next_unattenuated,
				   droplet_fraction[i1],
				   pristine_ice_fraction[i1],
				   afactor_near_AD[i1],
				   &width2g_AD,
				   &zeta2g_AD,
				   &covg_AD,
				   &Theta2_next_AD,
				   &bscat_next_unattenuated_AD,
				   &bscat_air_next_unattenuated_AD,
				   &(droplet_fraction_AD[i1]),
				   &(pristine_ice_fraction_AD[i1]));

	  ext_AD[i1] += bscat_next_unattenuated_AD/ext_bscat_ratio[i1];
	  // BUG FIX!!!
	  ext_bscat_ratio_AD[i1] -= bscat_next_unattenuated_AD
	    *bscat_next_unattenuated/ext_bscat_ratio[i1];
	  if (ext_air_AD) {
	    ext_air_AD[i1] += bscat_air_next_unattenuated_AD
	      *bscat_ext_ratio_air;
	  }
	  bscat_next_unattenuated_AD = 0.0;
	  bscat_air_next_unattenuated_AD = 0.0;

	  radius_AD[i1] -= Theta2_next_AD*Theta2_next*2.0/radius[i1];
	  Theta2_next_AD = 0.0;
	}
	// SECTION 2.3.1
	trans_AD[i1] += covg_AD*cov[i1]/transg;
	cov_AD[i1] += covg_AD*trans[i1]/transg;
	transu_AD[i1] -= covg_AD*covu/transg;
	transg_AD -= covg_AD*covg/transg;
	covg_AD = 0.0;

	trans_AD[i1] += zeta2g_AD*zeta2[i1]/transg;
	zeta2_AD[i1] += zeta2g_AD*trans[i1]/transg;
	transu_AD[i1] -= zeta2g_AD*rho_transmitter2/transg;
	transg_AD -= zeta2g_AD*zeta2g/transg;
	zeta2g_AD = 0.0;

	trans_AD[i1] += width2g_AD*width2[i1]/transg;
	width2_AD[i1] += width2g_AD*trans[i1]/transg;
	transu_AD[i1] -= width2g_AD*width2u/transg;
	transg_AD -= width2g_AD*width2g/transg;
	width2g_AD = 0.0;
	
	trans_AD[i1] += transg_AD;
	transu_AD[i1] -= transg_AD;
	transg_AD = 0.0;
      }
      // SECTION 2.2
      if (width2a_orig[i1] > width2a_max) {
	// SECTION 2.2.1
	/* The first Gaussian is larger than the receiver FOV and
	   needs adjusting */
	
	transa_AD[i1] += covr_AD[i1]*cova[i1]/transr[i1];
	cova_AD[i1] += covr_AD[i1]*transa[i1]/transr[i1];
	transu_AD[i1] += covr_AD[i1]*covu/transr[i1];
	transr_AD[i1] -= covr_AD[i1]*covr[i1]/transr[i1];
	covr_AD[i1] = 0.0;
     
	transa_AD[i1] += zeta2r_AD[i1]*zeta2a[i1]/transr[i1];
	zeta2a_AD[i1] += zeta2r_AD[i1]*transa[i1]/transr[i1];
	transu_AD[i1] += zeta2r_AD[i1]*rho_transmitter2/transr[i1];
	transr_AD[i1] -= zeta2r_AD[i1]*zeta2r[i1]/transr[i1];
	zeta2r_AD[i1] = 0.0;
	
	/*
	covr_orig_AD[i1] = covr_AD[i1];
	covr_AD[i1] = 0.0;

	zeta2r_orig_AD[i1] = zeta2r_AD[i1];
	zeta2r_AD[i1] = 0.0;

	width2r_orig_AD[i1] = width2r_AD[i1];
	width2r_AD[i1] = 0.0;
	*/
	
	transa_AD[i1] += width2r_AD[i1]*width2a[i1]/transr[i1];
	width2a_AD[i1] += width2r_AD[i1]*transa[i1]/transr[i1];
	transu_AD[i1] += width2r_AD[i1]*width2u/transr[i1];
	transr_AD[i1] -= width2r_AD[i1]*width2r[i1]/transr[i1];
	width2r_AD[i1] = 0.0;
	
	transu_AD[i1] += transr_AD[i1];
	transa_AD[i1] += transr_AD[i1];
	transr_AD[i1] = 0.0;

	ms_real correlation2a = cova_orig[i1]*cova_orig[i1]
	  /(width2a_orig[i1]*zeta2a_orig[i1]);
	ms_real factor = width2a_max/width2a_orig[i1];
	ms_real correlation2a_AD = zeta2a_AD[i1]*zeta2a_orig[i1]*(factor-1.0);
	ms_real factor_AD = zeta2a_AD[i1]*zeta2a_orig[i1]*correlation2a;
	zeta2a_orig_AD[i1] += zeta2a_AD[i1]
	  *(correlation2a*factor + 1.0-correlation2a);
	zeta2a_AD[i1] = 0.0;

	factor_AD += cova_AD[i1]*cova_orig[i1];
	cova_orig_AD[i1] += cova_AD[i1]*factor;
	cova_AD[i1] = 0.0;
	
	factor_AD += transa_AD[i1]*transa_orig[i1];
	transa_orig_AD[i1] += transa_AD[i1]*factor;
	transa_AD[i1] = 0.0;

	//transa_orig_AD[i1] = transa_AD[i1];
	//transa_AD[i1] = 0.0;

	width2a_AD[i1] = 0.0;
	//width2a_orig_AD[i1] = width2a_AD[i1];
	//width2a_AD[i1] = 0.0;
	width2a_orig_AD[i1] -= factor_AD*factor/width2a_orig[i1];
	factor_AD = 0.0;

	cova_orig_AD[i1] += correlation2a_AD*2.0*cova_orig[i1]
	  /(width2a_orig[i1]*zeta2a_orig[i1]);
	width2a_orig_AD[i1] -= correlation2a_AD*correlation2a/width2a_orig[i1];
	zeta2a_orig_AD[i1] -= correlation2a_AD*correlation2a/zeta2a_orig[i1];
	correlation2a_AD = 0.0;
      }
      else {
	// SECTION 2.2.2
	transa_orig_AD[i1] = transa_AD[i1];
	transa_AD[i1] = 0.0;

	width2a_orig_AD[i1] = width2a_AD[i1];
	width2a_AD[i1] = 0.0;

	cova_orig_AD[i1] = cova_AD[i1];
	cova_AD[i1] = 0.0;

	zeta2a_orig_AD[i1] = zeta2a_AD[i1];
	zeta2a_AD[i1] = 0.0;

	transr_orig_AD[i1] = transr_AD[i1];
	transr_AD[i1] = 0.0;

	width2r_orig_AD[i1] = width2r_AD[i1];
	width2r_AD[i1] = 0.0;

	covr_orig_AD[i1] = covr_AD[i1];
	covr_AD[i1] = 0.0;

	zeta2r_orig_AD[i1] = zeta2r_AD[i1];
	zeta2r_AD[i1] = 0.0;
      }

      // SECTION 2.1
      trans_AD[i1] += transb_AD[i1];
      transr_orig_AD[i1] -= transb_AD[i1];
      transb_AD[i1] = 0.0;

      transr_orig_AD[i1] += cova_orig_AD[i1]*covr_orig[i1]/transa_orig[i1];
      covr_orig_AD[i1] += cova_orig_AD[i1]*transr_orig[i1]/transa_orig[i1];
      transu_AD[i1] -= cova_orig_AD[i1]*covu/transa_orig[i1];
      transa_orig_AD[i1] -= cova_orig_AD[i1]*cova_orig[i1]/transa_orig[i1];
      cova_orig_AD[i1] = 0.0;

      transr_orig_AD[i1] += zeta2a_orig_AD[i1]*zeta2r_orig[i1]/transa_orig[i1];
      zeta2r_orig_AD[i1] += zeta2a_orig_AD[i1]*transr_orig[i1]/transa_orig[i1];
      transu_AD[i1] -= zeta2a_orig_AD[i1]*rho_transmitter2/transa_orig[i1];
      transa_orig_AD[i1] -= zeta2a_orig_AD[i1]*zeta2a_orig[i1]/transa_orig[i1];
      zeta2a_orig_AD[i1] = 0.0;
      
      transr_orig_AD[i1]
	+= width2a_orig_AD[i1]*width2r_orig[i1]/transa_orig[i1];
      width2r_orig_AD[i1]
	+= width2a_orig_AD[i1]*transr_orig[i1]/transa_orig[i1];
      transu_AD[i1] -= width2a_orig_AD[i1]*width2u/transa_orig[i1];
      transa_orig_AD[i1]
	-= width2a_orig_AD[i1]*width2a_orig[i1]/transa_orig[i1];
      width2a_orig_AD[i1] = 0.0;
    }
    else {
      // SECTION 5
      transr_orig_AD[i1] = transr_AD[i1];
      transr_AD[i1] = 0.0;
      
      width2r_orig_AD[i1] = width2r_AD[i1];
      width2r_AD[i1] = 0.0;
      
      covr_orig_AD[i1] = covr_AD[i1];
      covr_AD[i1] = 0.0;
      
      zeta2r_orig_AD[i1] = zeta2r_AD[i1];
      zeta2r_AD[i1] = 0.0;

      if (layer_od > 0.0) {
	// SECTION 3
	ms_real bscat_factor = transu[i]*(1.0-exp(-layer_od))/layer_od;
	ms_real bscat_factor_AD = 0.0;
	//	ms_real bscat_unattenuated_AD = 0.0;
	//	ms_real bscat_air_unattenuated_AD = 0.0;
	/* No forward scattering but some attenuation */
	if (bscat_air_out) {
	  bscat_factor_AD += bscat_AD[i]*bscat_unattenuated;
	  bscat_unattenuated_AD += bscat_AD[i]*bscat_factor;
	  bscat_factor_AD += bscat_air_AD[i]*bscat_air_unattenuated;
	  bscat_air_unattenuated_AD += bscat_air_AD[i]*bscat_factor;
	  //	bscat_AD[i] = 0.0;
	  //	bscat_air_AD[i] = 0.0;
	}
	else {
	  bscat_factor_AD += bscat_AD[i]
	    * (bscat_unattenuated+bscat_air_unattenuated);
	  bscat_unattenuated_AD += bscat_AD[i]*bscat_factor;
	  bscat_air_unattenuated_AD += bscat_AD[i]*bscat_factor;
	  //	bscat_AD[i] = 0.0;
	}
      
	layer_od_AD += bscat_factor_AD*(transu[i]/layer_od)
	  *(inv_exp_layer_od - (1.0-inv_exp_layer_od)/layer_od);
	transu_AD[i] += bscat_factor_AD*bscat_factor/transu[i];
	bscat_factor_AD = 0.0;
      }
      //    else {
      // SECTION 4
      //      bscat_AD[i] = 0.0;
      //      if (bscat_air_out) {
      //	bscat_air_AD[i] = 0.0;
      //      }
      //    }
    }
   
    // SECTION 1
    transr_orig_AD[i1] += transa_orig_AD[i1];
    transu_AD[i1] -= transa_orig_AD[i1];
    transa_orig_AD[i1] = 0.0;

    transr_AD[i] += transr_orig_AD[i1]*exp(-(ext[i]+2.0*ext_air[i])*drange);
    ext_AD[i] -= transr_orig_AD[i1]*transr_orig[i1]*drange;
    if (ext_air_AD) {
      ext_air_AD[i] -= transr_orig_AD[i1]*transr_orig[i1]*2.0*drange;
    }
    transr_orig_AD[i1] = 0.0;

    zeta2r_AD[i] += zeta2r_orig_AD[i1];
    ext_AD[i] += zeta2r_orig_AD[i1]*drange*Theta2;
    Theta2_AD += zeta2r_orig_AD[i1]*drange*ext[i];
    zeta2r_orig_AD[i1] = 0.0;

    covr_AD[i] += covr_orig_AD[i1];
    // FIXED drange2 -> drange
    zeta2r_AD[i] += covr_orig_AD[i1]*drange;
    // FIXED drange -> drange2
    ext_AD[i] += covr_orig_AD[i1]*Theta2*drange2*0.5;
    Theta2_AD += covr_orig_AD[i1]*ext[i]*drange2*0.5;
    covr_orig_AD[i1] = 0.0;

    width2r_AD[i] += width2r_orig_AD[i1];
    zeta2r_AD[i] += width2r_orig_AD[i1]*drange2;
    covr_AD[i] += width2r_orig_AD[i1]*2.0*drange;
    ext_AD[i] += width2r_orig_AD[i1]*Theta2*drange3*0.333333;
    Theta2_AD += width2r_orig_AD[i1]*ext[i]*drange3*0.333333;
    width2r_orig_AD[i1] = 0.0;
    
    trans_AD[i] += trans_AD[i1]*exp(-(ext[i]+2.0*ext_air[i])*drange);
    ext_AD[i] -= trans_AD[i1]*trans[i1]*drange;
    if (ext_air_AD) {
      ext_air_AD[i] -= trans_AD[i1]*trans[i1]*2.0*drange;
    }
    trans_AD[i1] = 0.0;

    zeta2_AD[i] += zeta2_AD[i1];
    ext_AD[i] += zeta2_AD[i1]*drange*Theta2;
    Theta2_AD += zeta2_AD[i1]*ext[i]*drange;
    zeta2_AD[i1] = 0.0;

    cov_AD[i] += cov_AD[i1];
    // FIXED drange2 -> drange
    zeta2_AD[i] += cov_AD[i1]*drange;
    // FIXED drange -> drange2
    ext_AD[i] += cov_AD[i1]*Theta2*drange2*0.5;
    Theta2_AD += cov_AD[i1]*ext[i]*drange2*0.5;
    cov_AD[i1] = 0.0;

    width2_AD[i] += width2_AD[i1];
    zeta2_AD[i] += width2_AD[i1]*drange2;
    cov_AD[i] += width2_AD[i1]*2.0*drange;
    ext_AD[i] += width2_AD[i1]*Theta2*drange3*0.333333;
    Theta2_AD += width2_AD[i1]*ext[i]*drange3*0.333333;
    width2_AD[i1] = 0.0;

    transu_AD[i] += transu_AD[i1]*inv_exp_layer_od;
    layer_od_AD -= transu_AD[i1]*transu[i1];
    transu_AD[i1] = 0.0;

    M_AD[i1] = 0.0;
    afactor_far_AD[i] = 0.0;
    afactor_near_AD[i1] = 0.0;

    // SECTION 0
    ext_AD[i] += layer_od_AD*2.0*drange;
    if (ext_air_AD) {
      ext_air_AD[i] += layer_od_AD*2.0*drange;
    }
    layer_od_AD = 0.0;

    if (ext_air_AD) {
      ext_air_AD[i] += bscat_air_unattenuated_AD*bscat_ext_ratio_air;
    }
    bscat_air_unattenuated_AD = 0.0;

    ext_AD[i] += bscat_unattenuated_AD/ext_bscat_ratio[i];
    ext_bscat_ratio_AD[i] -= bscat_unattenuated_AD
      *bscat_unattenuated/ext_bscat_ratio[i];
    bscat_unattenuated_AD = 0.0;
    
    radius_AD[i] -= Theta2_AD*2.0*Theta2/radius[i];
    Theta2_AD = 0.0;
  }


  return MS_SUCCESS;
}

