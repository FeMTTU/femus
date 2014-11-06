/* anisotropic_AD.c -- Adjoint of anisotropic factor for
   near-backscatter phase function

   Copyright (C) 2007-2010 Robin Hogan <r.j.hogan@reading.ac.uk> 

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

#define SQRT1_2 0.70710678118654752440
#define SQRT_PI 1.77245385090552

/* Numerical recipes function for the complementary error function */
static
ms_real
ms_erfc(ms_real x)
{
  ms_real t, z, ans;
  z = fabs(x);
  t = 1.0 / (1.0 + 0.5*z);
  ans = t * exp(-z*z - 1.26551223 
		+ t*(1.00002368 + 
		     t*(0.37409196 + 
			t*(0.09678418 + 
			   t*(-0.18628806 +
			      t*(0.27886807 +
				 t*(-1.13520398 +
				    t*(1.48851587 +
				       t*(-0.82215223 +
					  t*0.17087277)))))))));
  return x >= 0.0 ? ans : 2.0-ans;
}

/* Derivative of the complementary error function */
static 
inline
ms_real
ms_derfc_dx(ms_real x)
{
  return -exp(-x*x)*2.0/SQRT_PI;
}

/* Adjoint of the function to calculate the backscatter scaling factor
   due to an anisotropic near-backscatter phase function. The scaling
   factor is 1.0 for an isotropic phase function. The fraction of the
   cloud extinction (ext) that is due to ice and liquid in the volume
   may be specified separately; any remaining extinction is treated as
   isotropic. Theta2 is ignored in the ice parameterization, so should
   be set to the forward-lobe angular variance for droplets. Note that
   the near-backscatter phase function of ice particles is as
   predicted by Yang et al. (2000) whereas Baran (2004) asserts that
   an isotropic near-backscatter phase function is more realistic for
   irregular aggregates that dominate cirrus clouds; to implement this
   simply set ice_fraction to zero. */
ms_real
ms_anisotropic_factor_AD(
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
      ms_real* pristine_ice_fraction_AD)
{
  ms_real droplet_factor = 1.0;
  ms_real pristine_ice_factor = 1.0;
  ms_real factor;

  /* Calculate the variance of the scattering co-angle (pi minus the
     scattering angle) */
  ms_real orig_gamma2 = zeta2 + width2/(range*range) - 2.0*cov/range;
  ms_real gamma2 = orig_gamma2;
  /* If the variance of the distribution is larger than the receiver
     field-of-view then scale gamma2 appropriately */
  if (width2_max > 0.0 && width2 > width2_max) {
    ms_real cov_gamma = cov - width2/range;
    ms_real correlation2_gamma = cov_gamma*cov_gamma/(width2*zeta2);
    gamma2 = orig_gamma2
      * (correlation2_gamma*width2_max/width2 + 1.0 - correlation2_gamma);
  }

  /* Apply the 2-Gaussian parameterization for the Mie backscatter
     peak */
  if (droplet_fraction > 0.0) {
    droplet_factor = (1.0-MS_MIE_U1-MS_MIE_U2)
      + MS_MIE_U1/(1.0 + gamma2*MS_MIE_V1_SQD/Theta2)
      + MS_MIE_U2/(1.0 + gamma2*MS_MIE_V2_SQD/Theta2);
  }

  /* Apply the exponential parameterization for the Yang et al. ice
     particle phase functions */
  ms_real erfc = 0.0;
  ms_real sqrt_pi_exp_g2 = 0.0;
  ms_real sqrt_gamma2 = 0.0;
  if (pristine_ice_fraction > 0.0) {
    sqrt_gamma2 = sqrt(gamma2);
    ms_real g = sqrt_gamma2*0.5/MS_YANG_GAMMA0;
    sqrt_pi_exp_g2 = SQRT_PI*exp(g*g);
    erfc = ms_erfc(g);
    pristine_ice_factor = (1.0-MS_YANG_W)
      + MS_YANG_W*(1.0-g*sqrt_pi_exp_g2*erfc);
  }

  /* The final factor is a weighted average */
  ms_real weight
    = droplet_factor*droplet_fraction               /* Droplet contrib*/
      + pristine_ice_factor*pristine_ice_fraction   /* Ice contrib */
      + 1.0-droplet_fraction-pristine_ice_fraction; /* Isotropic remainder */

  factor
    = (bscat*weight 
       +bscat_air)                                  /* Isotropic air contrib */
    / (bscat + bscat_air);


  /* ADJOINT CALCULATION */

  *bscat_AD += factor_AD * (weight-factor) / (bscat+bscat_air);
  *bscat_air_AD += factor_AD * (1.0 - factor) / (bscat+bscat_air);
  ms_real weight_AD = factor_AD * bscat / (bscat+bscat_air);
  factor_AD = 0.0;

  ms_real droplet_factor_AD = weight_AD * droplet_fraction;
  *droplet_fraction_AD += weight_AD * (droplet_factor - 1.0);
  ms_real pristine_ice_factor_AD = weight_AD * pristine_ice_fraction;
  *pristine_ice_fraction_AD += weight_AD * (pristine_ice_factor - 1.0);
  weight_AD = 0.0;

  ms_real gamma2_AD = 0.0;
  if (pristine_ice_fraction > 0.0) {
    ms_real g = sqrt_gamma2*0.5/MS_YANG_GAMMA0;
    ms_real g_AD = -pristine_ice_factor_AD*MS_YANG_W*sqrt_pi_exp_g2
      *(erfc + 2.0*g*g*erfc + g*ms_derfc_dx(g));
    pristine_ice_factor_AD = 0.0;

    gamma2_AD += g_AD*0.25/(sqrt_gamma2*MS_YANG_GAMMA0);
    g_AD = 0.0;
  }

  if (droplet_fraction > 0.0) {
    ms_real denom1 = 1.0 + gamma2*MS_MIE_V1_SQD/Theta2;
    ms_real denom2 = 1.0 + gamma2*MS_MIE_V2_SQD/Theta2;
    gamma2_AD -= droplet_factor_AD
      *(MS_MIE_U1*MS_MIE_V1_SQD/(denom1*denom1*Theta2)
	+MS_MIE_U2*MS_MIE_V2_SQD/(denom2*denom2*Theta2));
    
    *Theta2_AD += droplet_factor_AD
      * (MS_MIE_U1*gamma2*MS_MIE_V1_SQD/(denom1*denom1*Theta2*Theta2)
	 +MS_MIE_U2*gamma2*MS_MIE_V2_SQD/(denom2*denom2*Theta2*Theta2));
    droplet_factor_AD = 0.0;
  }

  ms_real orig_gamma2_AD = 0.0;
  if (width2_max > 0.0 && width2 > width2_max) {
    ms_real cov_gamma = cov - width2/range;
    ms_real correlation2_gamma = cov_gamma*cov_gamma/(width2*zeta2);
    ms_real correlation2_gamma_AD = gamma2_AD*orig_gamma2
      * (width2_max/width2 - 1.0);
    *width2_AD -= gamma2_AD*orig_gamma2*correlation2_gamma
      /(width2*width2);
    orig_gamma2_AD += gamma2_AD 
      * (correlation2_gamma*width2_max/width2 + 1.0 - correlation2_gamma);
    gamma2_AD = 0.0;

    ms_real cov_gamma_AD = correlation2_gamma_AD*2.0*cov_gamma
      /(width2*zeta2);
    *width2_AD -= correlation2_gamma_AD*cov_gamma*cov_gamma
      /(width2*width2*zeta2);
    *zeta2_AD -= correlation2_gamma_AD*cov_gamma*cov_gamma
      /(width2*zeta2*zeta2);
    correlation2_gamma_AD = 0.0;

    *cov_AD += cov_gamma_AD;
    *width2_AD -= cov_gamma_AD/range;
    cov_gamma_AD = 0.0;
  }
  else {
    orig_gamma2_AD += gamma2_AD;
    gamma2_AD = 0.0;
  }
   
  *zeta2_AD += orig_gamma2_AD;
  *width2_AD += orig_gamma2_AD/(range*range);
  *cov_AD -= orig_gamma2_AD*2.0/range;
  orig_gamma2_AD = 0.0;

  return factor;
}

