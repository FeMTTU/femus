/* ms_wide_angle_AD.c -- Adjoint of wide angle calculation

   Copyright (C) 2010-2011 Robin Hogan <r.j.hogan@reading.ac.uk> 

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

/* Adjoint of the wide-angle part of the calculation */
int 
ms_wide_angle_AD(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,     /* Height of each range gate, metres */
    const ms_real *ext,       /* Total extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    const ms_real *ssa_air,
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    /* Adjoint inputs */
    const ms_real *bscat_AD,
    /* Adjoint outputs */
    ms_real *radius_AD,
    ms_real *ext_AD,
    ms_real *ssa_AD,
    ms_real *g_AD)
{
  ms_real transmittance[n]; /* Transmittance to each gate, including
			       internal attenuation */
  ms_real transmittance_mid[n+1]; /* Transmittance to each gate
				     edge */
  ms_real drange = fabs(range[2]-range[1]);
  int i, status;

  /* Local vectors allocated on the stack */
  /* Single scattering properties of clouds/aerosol and air, after
     delta-scaling if appropriate */
  ms_real ext_prime[n+1];
  ms_real ssa_prime[n+1];
  ms_real g_prime[n+1];
  /* Source terms in the inward and outward directions */
  ms_real src_power_in[n+1];
  ms_real src_power_out[n+1];
  /* Lateral variance of source photons (m2) */
  ms_real src_width2[n+1];

  /* Adjoint variables */
  ms_real ext_prime_AD[n+1];
  ms_real ssa_prime_AD[n+1];
  ms_real g_prime_AD[n+1];
  ms_real src_power_in_AD[n+1];
  ms_real src_power_out_AD[n+1];
  ms_real src_width2_AD[n+1];

  ms_real wide_factor[n];

  // FIX Insert check on field-of-view to optimize 
  /* Set each element of each vector to zero. */
  for (i = 0; i <= n; i++) {
    ext_prime[i] = 0.0;
    ssa_prime[i] = 1.0;
    g_prime[i] = 0.0;
    src_power_in[i] = 0.0;
    src_power_out[i] = 0.0;
    src_width2[i] = 0.0;
    ext_prime_AD[i] = ssa_prime_AD[i] = g_prime_AD[i] = 0.0;
    src_power_in_AD[i] = src_power_out_AD[i] = 0.0;
    src_width2_AD[i] = 0.0;
    wide_factor[i] = 1.0;
  }
  ext_prime_AD[n] = ssa_prime_AD[n] = g_prime_AD[n] = 0.0;
  src_power_in_AD[n] = src_power_out_AD[n] = 0.0;

  transmittance_mid[0] = 1.0;

  if (config->wide_angle_algorithm == MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE) {
    /* Section 1a */
    /* If there is no forward lobe then all scattered radiation should
       be considered by the wide-angle multiple scattering
       calculation */
    if (ext_air && ssa_air) {
      /* Section 1a.1 */
      /* Add the gas and particle contributions */
      for (i = 0; i < n; i++) {
	ext_prime[i] = ext_air[i] + ext[i];
	if (ext_prime[i] <= 0.0) {
	  /* Section 1a.1.1 */
	  ssa_prime[i] = 0.0;
	  g_prime[i] = 0.0;
	}
	else {
	  /* Section 1a.1.2 */
	  ssa_prime[i] = (ext_air[i]*ssa_air[i] + ext[i]*ssa[i])
	    / (ext_air[i]+ext[i]);
	  if (ssa_prime[i] <= 0.0) {
	    /* No scattering */
	    g_prime[i] = 0.0;
	  }
	  else {
	    /* The following assumes any molecular scattering to be
	       isotropic (i.e. its g is zero) */
	    g_prime[i] = ext[i]*ssa[i]*g[i]
	      /(ext_air[i]*ssa_air[i] + ext[i]*ssa[i]);
	  }
	}
      }
    }
    else {
      /* Section 1a.2 */
      for (i = 0; i < n; i++) {
	ext_prime[i] = ext[i];
	ssa_prime[i] = ssa[i];
	g_prime[i] = g[i];
      }
    }
    /* Calculate the source terms for the wide-angle calculation */
    for (i = 0; i < n; i++) {
      /* Section 1b */
      ms_real od = ext_prime[i]*drange; /* optical depth of layer */
      ms_real src_power;
      /* Section 1c */
      if (od > 0.0) {
	/* Section 1c.1 */
	transmittance_mid[i+1] = transmittance_mid[i]*exp(-od);
      }
      else {
	/* Section 1c.2 */
	transmittance_mid[i+1] 
	  = transmittance_mid[i];
      }
      /* Section 1d */
      /*
      src_power = config->coherent_backscatter_enhancement
	*config->ss_multiplier*transmittance_mid[i]
	*(1.0-exp(-ext_prime[i]*ssa_prime[i]*drange));
	*/
      src_power = config->coherent_backscatter_enhancement
	*config->ss_multiplier*ssa_prime[i]
	*(transmittance_mid[i]-transmittance_mid[i+1]);

      if (g_prime[i] >= 1.0/(3.0*MS_MU1)) {
	/* Section 1d.1 */
	/* If g exceeds 2/3, all energy is sent into outward diffuse stream */
	src_power_in[i] = 0.0;
	src_power_out[i] = src_power;
      }
      else if (g_prime[i] <= -1.0/(3.0*MS_MU1)) {
	/* Section 1d.2 */
	/* If g is less than -2/3, all energy is sent into inward
	   diffuse stream (this never happens in atmospheric
	   scattering) */
	src_power_in[i] = src_power;
	src_power_out[i] = 0.0;
      }
      else {
	/* Section 1d.3 */
	/* Two-stream phase function */
	src_power_in[i] = src_power*0.5*(1.0-3.0*g_prime[i]*MS_MU1);
	src_power_out[i] = src_power*0.5*(1.0+3.0*g_prime[i]*MS_MU1);
      }
      /* Section 1e */
      src_width2[i] = instrument.rho_transmitter*instrument.rho_transmitter
	*(range[i]-instrument.altitude)*(range[i]-instrument.altitude);
    } /* End of spatial loop */
  }

  else {
    /* Section 2a */
    /* Calculate the source terms for the wide-angle calculation in
       the case that there is a forward lobe and so some of the
       scattered radiation is not treated by the wide-angle part of
       the code */
    
    if (config->options & MS_IGNORE_SOURCE_WIDENING) {
      /* Include only the contribution of forward scattering to the
	 power, not to the width */
      for (i = 0; i < n; i++) {
	src_width2[i] = instrument.rho_transmitter*instrument.rho_transmitter
	  *(range[i]-instrument.altitude)*(range[i]-instrument.altitude);
      }
    }
    else {
      /* Calculate the contribution to the width */
      ms_variance(n, instrument.wavelength, instrument.rho_transmitter,
		  instrument.altitude,
		  range, radius, ext, src_width2);
    }

    /* Section 2b */
    for (i = 0; i < n; i++) {
      ms_real od, od_factor, power_factor, src_power, g_equiv;

      if (config->options & MS_SSA_SCALES_FORWARD_LOBE) {
	/* Section 2b.1 */
	/* Single-scattering albedo of less than one acts to reduce the
	   phase function for all angles; this means that both the
	   forward lobe enhancement to the transmission is reduced, and
	   the side-scattered power is reduced */
	/* Note that currently this option is not treated at all in the
	   ms_small_angle() function */
	od_factor = 2.0-ssa[i];
	power_factor = ssa[i];
      }
      else if (ssa[i] >= 0.5/g[i] || ext[i] <= 0.0) {
	/* Section 2b.2 */
	/* For single-scattering albedo of greater than 0.5/g, assume all
	   the absorption affects only the wide-angle part of the phase
	   function, not the forward lobe */
	/* If the extinction is zero then ssa, g, od_factor and
	   power_factor don't matter so we may as well do the same
	   thing here even if ssa=g=0 */
	od_factor = 1.0;
	power_factor = 2.0*ssa[i]-1.0;
      }
      else {
	/* Section 2b.3 */
	/* For single-scattering albedo less than 0.5/g,
	   "diffraction-scaled" value of g becomes negative - this
	   doesn't make sense */
	fprintf(stderr, "Error: cannot simulate multiple scattering including a forward lobe when ssa < 0.5/g and ext > 0 (ssa=%g, g=%g)\n", ssa[i], g[i]);
	return MS_FAILURE;
	//	od_factor = 2.0*(1.0-ssa[i]);
	//	power_factor = 0.0;
      }

      /* Section 2c */
      /* Calculate the optical depth to the start of gate i */
      if (ext_air) {
	/* Section 2c.1 */
	od = (0.5*ext[i]*od_factor + ext_air[i])*drange;
      }
      else {
	/* Section 2c.2 */
	od = 0.5*ext[i]*od_factor*drange;
      }
      /* Section 2d */
      /* Calculate the transmittance to gate i, including in-gate
	 attenuation */
      if (od > 0.0) {
	/* Section 2d.1a */
	transmittance_mid[i+1] = transmittance_mid[i]*exp(-od);
	transmittance[i] = (transmittance_mid[i]-transmittance_mid[i+1])
	  / od;
	/* Section 2d.1b */
	if (ext_air && ssa_air) {
	  /* Section 2d.1b.1 */
	  wide_factor[i] = (ext[i]*ssa[i] + ext_air[i]*ssa_air[i])
	    /(0.5*ext[i]*power_factor + ext_air[i]*ssa_air[i]);
	}
	else {
	  /* Section 2d.1b.2 */
	  wide_factor[i] = 2.0*ssa[i]/power_factor;
	}
      }
      else {
	/* Section 2d.2 */
	transmittance_mid[i+1]
	  = transmittance[i] = transmittance_mid[i];
	wide_factor[i] = 1.0;
      }
      
      /* Section 2e */
      /* Calculate source power and source width-squared for the
	 wide-angle calculation */
      if (ext_air && ssa_air) {
	/* Section 2e.1 */
	src_power = config->coherent_backscatter_enhancement
	  *config->ss_multiplier*transmittance[i]
	  *(0.5*ext[i]*power_factor + ext_air[i]*ssa_air[i])*drange;
      }
      else {
	/* Section 2e.2 */
	src_power = config->coherent_backscatter_enhancement
	  *config->ss_multiplier*transmittance[i]
	  *0.5*ext[i]*power_factor*drange;
      }
      /* Section 2f */
      g_equiv = 1.0-wide_factor[i]+wide_factor[i]*g[i];
      /* Section 2g */
      if (g_equiv >= 1.0/(3.0*MS_MU1)) {
	/* Section 2g.1 */
	src_power_in[i] = 0.0;
	src_power_out[i] = src_power;
      }
      else {
	/* Section 2g.2 */
	src_power_in[i] = src_power*0.5*(1.0-3.0*g_equiv*MS_MU1);
	src_power_out[i] = src_power*0.5*(1.0+3.0*g_equiv*MS_MU1);
      }
    } /* End of spatial loop */
    
    /* Section 2h */
    if (ext_air && ssa_air) {
      /* Weighted delta scaling assuming that the gaseous asymmetry
	 factor is 0.0 */
      for (i = 0; i < n; i++) {
	/* Section 2h.1a */
	ext_prime[i] = ext_air[i] + ext[i]*(1.0-ssa[i]*g[i]*g[i]);
	/* Section 2h.1b */
	if (ext_air[i] + ext[i] <= 0.0) {
	  /* Section 2h.1b.1 */
	  ssa_prime[i] = 0.0;
	  g_prime[i] = 0.0;
	}
	else {
	  /* Section 2h.1b.2 */
	  ssa_prime[i] 
	    = (ext_air[i]*ssa_air[i]
	       + ext[i]*ssa[i]*(1.0-g[i]*g[i]))
	    / (ext_air[i]+ext[i]*(1.0-ssa[i]*g[i]*g[i]));
	  /* Assume g of air is zero (isotropic scattering) */
	  g_prime[i] = ext[i]*ssa[i]*g[i]/(1.0+g[i])
	    / (ext_air[i]*ssa_air[i]/(1.0-g[i]*g[i])
	       + ext[i]*ssa[i]);
	}
      }
    }
    else {
      /* Section 2h.2 */
      /* Standard delta scaling without air contribution */
      for (i = 0; i < n; i++) {
	ext_prime[i] = ext[i]*(1.0-ssa[i]*g[i]*g[i]);
	ssa_prime[i] = ssa[i]*(1.0-g[i]*g[i])/(1.0-ssa[i]*g[i]*g[i]);
	g_prime[i] = g[i]/(1.0+g[i]);
      }
    }
  }

  /* Perform the 2-stream wide-angle multiple scattering
     calculation, including adjoint */

  /*
  fprintf(stderr, "ext'  ssa'   g'   src_power_in   src_power_out   src_width2\n");
  for (i = 0; i < n; i++) {
    fprintf(stderr, "%10g %10g %10g %10g %10g %10g\n", ext_prime[i],
	    ssa_prime[i], g_prime[i], src_power_in[i], src_power_out[i],
	    src_width2[i]);
  }
  */  

  status = 
    ms_tdts_AD(n, m, config, instrument, surface, 
	       range, ext_prime, ssa_prime, g_prime,
	       src_power_in, src_power_out, src_width2,
	       bscat_out,
	       bscat_AD,
	       ext_prime_AD, ssa_prime_AD, g_prime_AD,
	       src_power_in_AD, src_power_out_AD,
	       src_width2_AD);
  /*
  fprintf(stderr, "ext'_AD  ssa'_AD   g'_AD   src_power_in_AD   src_power_out_AD   src_width2_AD\n");
  for (i = 0; i < n; i++) {
    fprintf(stderr, "%10g %10g %10g %10g %10g %10g\n", ext_prime_AD[i],
	    ssa_prime_AD[i], g_prime_AD[i], src_power_in_AD[i], src_power_out_AD[i],
	    src_width2_AD[i]);
  }
  */

  if (status != MS_SUCCESS) {
    return status;
  }

  /* Main adjoint code */
  if (config->wide_angle_algorithm == MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE) {
    /* If there is no forward lobe then all scattered radiation should
       be considered by the wide-angle multiple scattering
       calculation */
    ms_real transmittance_mid_next_AD = 0.0;
    ms_real transmittance_mid_AD = 0.0;

    /* Calculate the source terms for the wide-angle calculation */
    for (i = n-1; i >= 0; i--) {
      ms_real od = ext_prime[i]*drange; /* optical depth of layer */
      ms_real src_power_AD, od_AD = 0.0;
      /* No section 1e: we don't need to do adjoint of src_width2 */
      if (g_prime[i] >= 1.0/(3.0*MS_MU1)) {
	/* Section 1d.1 */
	src_power_AD = src_power_out_AD[i];
      }
      else if (g_prime[i] <= -1.0/(3.0*MS_MU1)) {
	/* Section 1d.2 */
	/* If g is less than -2/3, all energy is sent into inward
	   diffuse stream (this never happens in atmospheric
	   scattering) */
	src_power_AD = src_power_in_AD[i];
      }
      else {
	/* Section 1d.3 */
	ms_real src_power = src_power_in[i] + src_power_out[i];
	/* Two-stream phase function */
	src_power_AD = 0.5*((src_power_in_AD[i]*(1.0-3.0*g_prime[i]*MS_MU1))
			    +(src_power_out_AD[i]*(1.0+3.0*g_prime[i]*MS_MU1)));
	g_prime_AD[i] += src_power*1.5*MS_MU1
	  *(src_power_out_AD[i]-src_power_in_AD[i]);
      }
      src_power_out_AD[i] = src_power_in_AD[i] = 0.0;

      /* Section 1d */
      /*
      ms_real exp_factor = exp(-ext_prime[i]*ssa_prime[i]*drange);
      transmittance_mid_AD += config->coherent_backscatter_enhancement
	*config->ss_multiplier*(1.0-exp_factor)*src_power_AD;
      ext_prime_AD[i] += config->coherent_backscatter_enhancement
	*config->ss_multiplier*transmittance_mid[i]*exp_factor
	*ssa_prime[i]*drange*src_power_AD;
      ssa_prime_AD[i] += config->coherent_backscatter_enhancement
	*config->ss_multiplier*transmittance_mid[i]*exp_factor
	*ext_prime[i]*drange*src_power_AD;
      src_power_AD = 0.0;
      */
      ssa_prime_AD[i] += src_power_AD
	*config->coherent_backscatter_enhancement*config->ss_multiplier
	*(transmittance_mid[i]-transmittance_mid[i+1]);
      transmittance_mid_AD += src_power_AD
	*config->coherent_backscatter_enhancement*config->ss_multiplier
	*ssa_prime[i];
      transmittance_mid_next_AD -= src_power_AD
	*config->coherent_backscatter_enhancement*config->ss_multiplier
	*ssa_prime[i];
      src_power_AD = 0.0;      

      /* Section 1c */
      if (od > 0.0) {
	/* Section 1c.1 */
	transmittance_mid_AD += exp(-od)*transmittance_mid_next_AD;
	od_AD -= transmittance_mid[i+1]*transmittance_mid_next_AD;
      }
      else {
	/* Section 1c.2 */
	transmittance_mid_AD += transmittance_mid_next_AD;
      }
      transmittance_mid_next_AD = 0.0;

      /* Section 1b */
      ext_prime_AD[i] += drange*od_AD;
      od_AD = 0.0;

      /* Redefine the adjoints of the mid-point transmittances so that
	 we don't need to keep an array of values */
      transmittance_mid_next_AD = transmittance_mid_AD;
      transmittance_mid_AD = 0.0;
    } /* End of spatial loop */

    if (ext_air && ssa_air) {
      /* Section 1a.1 */
      /* Add the gas and particle contributions */
      for (i = 0; i < n; i++) {
	if (ext_prime[i] <= 0.0) {
	  /* Section 1a.1.1 */
	  ssa_prime_AD[i] = 0.0;
	  g_prime_AD[i] = 0.0;
	}
	else {
	  /* Section 1a.1.2 */
	  if (ssa_prime[i] <= 0.0) {
	    g_prime_AD[i] = 0.0;
	  }
	  else {
	    /* The following assumes any molecular scattering to be
	       isotropic (i.e. its g is zero) */
	    ms_real factor = 1.0/(ext_air[i]*ssa_air[i] + ext[i]*ssa[i]);
	    g_AD[i] += factor*ext[i]*ssa[i]*g_prime_AD[i];
	    ext_AD[i] += factor*ssa[i]*g[i]*g_prime_AD[i]
	      *(1.0 - ext[i]*ssa[i]*factor);
	    ssa_AD[i] += factor*ext[i]*g[i]*g_prime_AD[i]
	      *(1.0 - ext[i]*ssa[i]*factor);
	    g_prime_AD[i] = 0.0;
	  }
	  ms_real factor_ext = 1.0/(ext_air[i] + ext[i]);
	  ext_AD[i] += factor_ext*ssa_prime_AD[i]
	    *(ssa[i] - ssa_prime[i]);
	  ssa_AD[i] += factor_ext*ext[i]*ssa_prime_AD[i];
	  ssa_prime_AD[i] = 0.0;
	}
	ext_AD[i] += ext_prime_AD[i];
	ext_prime_AD[i] = 0.0;
      }
    }
    else {
      /* Section 1a.2 */
      for (i = 0; i < n; i++) {
	ext_AD[i] += ext_prime_AD[i];
	ssa_AD[i] += ssa_prime_AD[i];
	g_AD[i] += g_prime_AD[i];
	ext_prime_AD[i] = 0.0;
	ssa_prime_AD[i] = 0.0;
	g_prime_AD[i] = 0.0;
      }
    }
  }

  else {
    /* Calculate the source terms for the wide-angle calculation in
       the case that there is a forward lobe and so some of the
       scattered radiation is not treated by the wide-angle part of
       the code */
    ms_real transmittance_mid_next_AD = 0.0;
    ms_real transmittance_mid_AD = 0.0;

    /* Section 2h */
    if (ext_air && ssa_air) {
      // BUG IN THIS BLOCK? 
      /* Weighted delta scaling assuming that the gaseous asymmetry
	 factor is 0.0 */
      for (i = 0; i < n; i++) {
	/* Section 2h.1b */
	if (ext_air[i] + ext[i] <= 0.0) {
	  /* Section 2h.1b.1 */
	  ssa_prime_AD[i] = 0.0;
	  g_prime_AD[i] = 0.0;
	}
	else {
	  /* Section 2h.1b.2 */
	  ms_real factor = 1.0/(ext_air[i]*ssa_air[i]/(1.0-g[i]*g[i])
				+ ext[i]*ssa[i]);
	  ms_real factor_ssa = 1.0/(ext_air[i]
				    + ext[i]*(1.0-ssa[i]*g[i]*g[i]));
	  g_AD[i] += g_prime[i]*g_prime_AD[i]
	    *(1.0/g[i] - 1.0/(1.0+g[i])
	      - 2.0*g[i]*ext_air[i]*ssa_air[i]*factor
	      /((1.0-g[i]*g[i])*(1.0-g[i]*g[i])));
	  ssa_AD[i] += factor*ext[i]*g_prime_AD[i]*g[i]/(1.0+g[i])
	    *(1.0 - factor*ext[i]*ssa[i]);
	  ext_AD[i] += factor*ssa[i]*g_prime_AD[i]*g[i]/(1.0+g[i])
	    *(1.0 - factor*ext[i]*ssa[i]);
	  g_prime_AD[i] = 0.0;

	  ext_AD[i] += factor_ssa*ssa_prime_AD[i]
	    *(ssa[i]*(1.0-g[i]*g[i]) - ssa_prime[i]*(1.0-ssa[i]*g[i]*g[i]));
	  ssa_AD[i] += factor_ssa*ssa_prime_AD[i]*ext[i]
	    *(1.0 + (ssa_prime[i]-1.0)*g[i]*g[i]);
	  g_AD[i] += factor_ssa*ssa_prime_AD[i]*2.0*g[i]*ext[i]*ssa[i]
	    *(ssa_prime[i] - 1.0);
	  ssa_prime_AD[i] = 0.0;
	}

	/* Section 2h.1a */
	ext_AD[i] += (1.0-ssa[i]*g[i]*g[i])*ext_prime_AD[i];
	ssa_AD[i] -= ext[i]*g[i]*g[i]*ext_prime_AD[i];
	g_AD[i] -= 2.0*ext[i]*ssa[i]*g[i]*ext_prime_AD[i];
	ext_prime_AD[i] = 0.0;
      }

    }
    else {
      /* Section 2h.2 */
      /* Standard delta scaling without air contribution */
      for (i = 0; i < n; i++) {
	ms_real factor_g = 1.0/(1.0+g[i]);
	ms_real factor_ssa = 1.0/(1.0-ssa[i]*g[i]*g[i]);
	g_AD[i] += factor_g*(1.0 - g[i]*factor_g)*g_prime_AD[i];
	g_prime_AD[i] = 0.0;
	ssa_AD[i] += (1.0-g[i]*g[i])*factor_ssa
	  *(1.0 + g[i]*g[i]*ssa[i]*factor_ssa)*ssa_prime_AD[i];
	g_AD[i] += ssa[i]*factor_ssa*2.0*g[i]
	  *((1.0-g[i]*g[i])*factor_ssa*ssa[i] - 1.0)*ssa_prime_AD[i];
	ssa_prime_AD[i] = 0.0;
	ext_AD[i] += (1.0-ssa[i]*g[i]*g[i])*ext_prime_AD[i];
	ssa_AD[i] -= ext[i]*g[i]*g[i]*ext_prime_AD[i];
	g_AD[i] -= 2.0*ext[i]*ssa[i]*g[i]*ext_prime_AD[i];
	ext_prime_AD[i] = 0.0;
      }
    }

    for (i = n-1; i >= 0; i--) {
      ms_real od;
      ms_real od_factor, power_factor;
      ms_real src_power = src_power_in[i]+src_power_out[i];
      ms_real g_equiv = 1.0-wide_factor[i]+wide_factor[i]*g[i];
      ms_real g_equiv_AD = 0.0, wide_factor_AD = 0.0, od_AD = 0.0;
      ms_real od_factor_AD = 0.0;
      ms_real transmittance_AD = 0.0;
      ms_real src_power_AD = 0.0;
      ms_real power_factor_AD = 0.0;

      /* Repeat sections 2b and 2c (not the adjoint part) to get od_factor,
	 power_factor and od */
      if (config->options & MS_SSA_SCALES_FORWARD_LOBE) {
	/* Section 2b.1 */
	od_factor = 2.0-ssa[i];
	power_factor = ssa[i];
      }
      else if (ssa[i] >= 0.5/g[i] || ext[i] <= 0.0) {
	/* Section 2b.2 */
	od_factor = 1.0;
	power_factor = 2.0*ssa[i]-1.0;
      }
      else {
	// Should not get here!
	fprintf(stderr, "Error: cannot simulate multiple scattering including a forward lobe when ssa < 0.5/g and ext > 0 (ssa=%g, g=%g)\n", ssa[i], g[i]);
	return MS_FAILURE;
	/* Section 2b.3 */
	//	od_factor = 2.0*(1.0-ssa[i]);
	//	power_factor = 0.0;
      }
      /* Section 2c */
      /* Calculate the optical depth to the start of gate i */
      if (ext_air) {
	/* Section 2c.1 */
	od = (0.5*ext[i]*od_factor + ext_air[i])*drange;
      }
      else {
	/* Section 2c.2 */
	od = 0.5*ext[i]*od_factor*drange;
      }

      /* Section 2g */
      if (g_equiv >= 1.0/(3.0*MS_MU1)) {
	/* Section 2g.1 */
	src_power_AD = src_power_out_AD[i];
      }
      else {
	/* Section 2g.2 */
	src_power_AD += 0.5*src_power_out_AD[i]*(1.0+3.0*g_equiv*MS_MU1)
	  + 0.5*src_power_in_AD[i]*(1.0-3.0*g_equiv*MS_MU1);
	g_equiv_AD = 1.5*src_power*MS_MU1
	  *(src_power_out_AD[i] - src_power_in_AD[i]);
      }
      src_power_out_AD[i] = src_power_in_AD[i] = 0.0;
      /* Section 2f */
      g_AD[i] += wide_factor[i]*g_equiv_AD;
      wide_factor_AD += (g[i]-1.0)*g_equiv_AD;
      g_equiv_AD = 0.0;

      /* Calculate source power and source width-squared for the
	 wide-angle calculation */
      /* Section 2e */
      if (ext_air && ssa_air) {
	/* Section 2e.1 */
	transmittance_AD += config->coherent_backscatter_enhancement
	  *config->ss_multiplier*src_power_AD
	  *(0.5*ext[i]*power_factor + ext_air[i]*ssa_air[i])*drange;
	ext_AD[i] += 0.5*config->coherent_backscatter_enhancement
	  *config->ss_multiplier*transmittance[i]
	  *power_factor*drange*src_power_AD;
	power_factor_AD += 0.5*config->coherent_backscatter_enhancement
	  *config->ss_multiplier*transmittance[i]
	  *ext[i]*drange*src_power_AD;
      }
      else {
	/* Section 2e.2 */
	ms_real factor = 0.5*config->coherent_backscatter_enhancement
	  *config->ss_multiplier*drange*src_power_AD;
	transmittance_AD += ext[i]*power_factor*factor;
	ext_AD[i] += transmittance[i]*power_factor*factor;
	power_factor_AD += ext[i]*transmittance[i]*factor;
      }
      src_power_AD = 0.0;
      
      if (od > 0.0) {
	if (ext_air && ssa_air) {
	  /* Section 2d.1b.1 */
	  ms_real factor = 1.0/(0.5*ext[i]*power_factor 
				+ ext_air[i]*ssa_air[i]);
	  ext_AD[i] += wide_factor_AD*factor
	    *(ssa[i] - 0.5*power_factor*wide_factor[i]);
	  ssa_AD[i] += ext[i]*factor*wide_factor_AD;
	  power_factor_AD -= factor*wide_factor[i]
	    *0.5*ext[i]*wide_factor_AD;
	}
	else {
	  /* Section 2d.1b.2 */
	  ssa_AD[i] += 2.0*wide_factor_AD/power_factor;
	  power_factor_AD -= 2.0*ssa[i]*wide_factor_AD
	    /(power_factor*power_factor);
	}
	wide_factor_AD = 0.0;

	/* Section 2d.1a */
	transmittance_mid_AD += transmittance_AD/od;
	transmittance_mid_next_AD -= transmittance_AD/od;
	od_AD -= transmittance[i]*transmittance_AD/od;
	transmittance_AD = 0.0;
	transmittance_mid_AD += exp(-od)*transmittance_mid_next_AD;
	od_AD -= transmittance_mid[i+1]*transmittance_mid_next_AD;
      }
      else {
	/* Section 2d.2 */
	wide_factor_AD = 0.0;
	/* CHECK Should this be += ? */
	transmittance_mid_AD = transmittance_mid_next_AD;
	transmittance_AD = transmittance_mid_next_AD;
      }
      transmittance_mid_next_AD = 0.0;

      /* Section 2c (merge of original sections 2c.1 and 2c.2 */
      ext_AD[i] += 0.5*od_factor*drange*od_AD;
      od_factor_AD += 0.5*ext[i]*drange*od_AD;
      od_AD = 0.0;

      if (config->options & MS_SSA_SCALES_FORWARD_LOBE) {
	/* Section 2b.1 */
	/* Single-scattering albedo of less than one acts to reduce the
	   phase function for all angles; this means that both the
	   forward lobe enhancement to the transmission is reduced, and
	   the side-scattered power is reduced */
	/* Note that currently this option is not treated at all in the
	   ms_small_angle() function */
	ssa_AD[i] += power_factor_AD - od_factor_AD;
	power_factor_AD = 0.0;
	od_factor_AD = 0.0;
      }
      else if (ssa[i] > 0.5) {
	/* Section 2b.2 */
	/* For single-scattering albedo of greater than 0.5, assume all
	   the absorption affects only the wide-angle part of the phase
	   function, not the forward lobe */
	ssa_AD[i] += 2.0*power_factor_AD;
	power_factor_AD = 0.0;
	od_factor_AD = 0.0;
      }
      else {
	/* Section 2b.3 */
	/* For single-scattering albedo less than 0.5, assume that there
	   is no wide-angle scattering and now the forward lobe starts
	   to be reduced */
	power_factor_AD = 0.0;
	ssa_AD[i] -= 2.0*od_factor_AD;
	od_factor_AD = 0.0;
      }
      /* Redefine the adjoints of the mid-point transmittances so that
	 we don't need to keep an array of values */
      transmittance_mid_next_AD = transmittance_mid_AD;
      transmittance_mid_AD = 0.0;
    } /* End of spatial loop */

    /* Section 2a */
    if (!(config->options & MS_IGNORE_SOURCE_WIDENING)) {
      ms_variance_AD(n, instrument.wavelength, instrument.rho_transmitter,
		     instrument.altitude,
		     range, radius, ext, src_width2, src_width2_AD,
		     radius_AD, ext_AD);
    }
  }

  return MS_SUCCESS;
}
