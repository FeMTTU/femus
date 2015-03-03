/* anisotropic.c -- Account for anisotropic near-backscatter phase function

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

#define SQRT1_2 0.70710678118654752440
#define SQRT_PI 1.77245385090552

/* Numerical recipes function for the complementary error function */
static
areal
ms_erfc(areal x)
{
  areal t, z, ans;
  z = abs(x);
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
  //  return x >= 0.0 ? ans : 2.0-ans;
  if (x >= 0.0) {
    return ans;
  }
  else {
    return 2.0-ans;
  }
}

/* Calculate the backscatter scaling factor due to an anisotropic
   near-backscatter phase function. This is 1.0 for an isotropic phase
   function. The fraction of the cloud extinction (ext) that is due to
   ice and liquid in the volume may be specified separately; any
   remaining extinction is treated as isotropic. Theta2 is ignored in
   the ice parameterization, so should be set to the forward-lobe
   angular variance for droplets. Note that the near-backscatter phase
   function of ice particles is as predicted by Yang et al. (2000)
   whereas Baran (2004) asserts that an isotropic near-backscatter
   phase function is more realistic for irregular aggregates that
   dominate cirrus clouds; to implement this simply set ice_fraction
   to zero. */
areal
ms_anisotropic_factor(
      areal width2,          /* spatial variance (m2) */
      areal zeta2,           /* angular variance (rad2) */
      areal cov,             /* covariance (m rad) */
      areal Theta2,          /* forward lobe variance (rad2) */
      ms_real range,           /* range from instrument (m) */
      ms_real width2_max,      /* receiver FOV variance (m2) */
      areal bscat,           /* Cloud bscat coefficient (m-1) */
      ms_real bscat_air,       /* Air bscat coefficient (m-1) */
      areal droplet_fraction,/* Fraction of ext from droplets */
      areal pristine_ice_fraction /* Fraction of ext from pristine ice */
		      )
{
  areal droplet_factor = 1.0;
  areal pristine_ice_factor = 1.0;
  areal factor;

  /* Calculate the variance of the scattering co-angle (pi minus the
     scattering angle) */
  areal gamma2 = zeta2 + width2/(range*range) - 2.0*cov/range;

  /* If the variance of the distribution is larger than the receiver
     field-of-view then scale gamma2 appropriately */
  if (width2_max > 0.0 && width2 > width2_max) {
    areal cov_gamma = cov - width2/range;
    areal correlation2_gamma = cov_gamma*cov_gamma/(width2*zeta2);
    gamma2 = gamma2
      *(correlation2_gamma*width2_max/width2 + 1.0 - correlation2_gamma);
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
  if (pristine_ice_fraction > 0.0) {
    areal g = sqrt(gamma2)*0.5/MS_YANG_GAMMA0;
    pristine_ice_factor = (1.0-MS_YANG_W)
      + MS_YANG_W*(1.0-g*SQRT_PI*exp(g*g)*ms_erfc(g));
  }

  /* The final factor is a weighted average */
  factor
    = (bscat*(droplet_factor*droplet_fraction      /* Droplet contrib*/
	      + pristine_ice_factor*pristine_ice_fraction /* Ice contrib */
	      + 1.0-droplet_fraction-pristine_ice_fraction)/* Isotropic
							      remainder */
       +bscat_air                                  /* Isotropic air contrib */
       ) / (bscat+bscat_air);

  return factor;
}

