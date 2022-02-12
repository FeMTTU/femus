/* wide_angle.c -- Wide-angle multiple scattering calculation

   Copyright (C) 2006-2011 Robin Hogan <r.j.hogan@reading.ac.uk> 

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

/* Do the wide-angle part of the calculation */
int ms_wide_angle(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const areal *radius,     /* Height of each range gate, metres */
    const areal *ext,       /* Total extinction coefficient, m-1 */
    const areal *ssa,       /* Total single-scatter albedo */
    const areal *g,         /* Total asymmetry factor */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    const ms_real *ssa_air,
    /* Output data */
    areal *bscat_out)       /* Measured backscatter, m-1 sr-1 */
{
  areal optical_depth = 0.0;
  ms_real drange = fabs(range[2]-range[1]);
  int i;

  /* Local vectors allocated on the stack */
  /* Single scattering properties of clouds/aerosol and air, after
     delta-scaling if appropriate */
  areal ext_prime[n+1];
  areal ssa_prime[n+1];
  areal g_prime[n+1];
  /* Source terms in the inward and outward directions */
  areal src_power_in[n+1];
  areal src_power_out[n+1];
  /* Lateral variance of source photons (m2) */
  areal src_width2[n+1];

/*   int first_gate = config->first_wide_angle_gate; */
/*   if (first_gate < 0) { */
/*     /\* Compare mean-free-path to footprint to work out when to start */
/*        wide-angle calculation *\/ */
/*     for (i = 0; i < n-1; i++) { */
/*       if (ext[i] > 0.0 && ssa[i] > 0.0) { */
/* 	/\* The footprint radius at the altitude of the cloud of the */
/* 	   widest field of view *\/ */
/* 	ms_real footprint = fabs(instrument.altitude-range[i]) */
/* 	  *instrument.rho_receiver[instrument.nfov-1]; */
/* 	/\* The relevant mean-free-path *\/ */
/* 	ms_real mfp = 1.0/(ext[i]*ssa[i]*(1.0 - g[i])); */
/* 	if (mfp < footprint*MS_MFP_FOOTPRINT_THRESHOLD) { */
/* 	  /\* Threshold met: start multiple scattering calculation from */
/* 	     this gate *\/ */
/* 	  first_gate = i; */
/* 	  break; */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   if (first_gate < 0) { */
/*     /\* No gate with small enough mean-free-path found: skip wide-angle */
/*        calculation all together *\/ */
/*     return MS_SUCCESS; */
/*   } */

  /* Set each element of each vector to zero. */
  for (i = 0; i <= n; i++) {
    ext_prime[i] = 0.0;
    ssa_prime[i] = 1.0;
    g_prime[i] = 0.0;
    src_power_in[i] = 0.0;
    src_power_out[i] = 0.0;
    src_width2[i] = 0.0;

#ifdef CPPAD
    /*
    ms_real tmp[4];
    tmp[0] = Value(Var2Par(radius[i]));
    tmp[1] = Value(Var2Par(ext[i]));
    tmp[2] = Value(Var2Par(ssa[i]));
    tmp[3] = Value(Var2Par(g[i]));
    fprintf(stderr, "??? %d %g %g %g %g\n", i, tmp[0], tmp[1], tmp[2], tmp[3]);
    */
#endif

  }

  if (config->wide_angle_algorithm == MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE) {
    /* If there is no forward lobe then all scattered radiation should
       be considered by the wide-angle multiple scattering
       calculation */
    if (ext_air && ssa_air) {
      /* Add the gas and particle contributions */
      for (i = 0; i < n; i++) {
	ext_prime[i] = ext_air[i] + ext[i];
	if (ext_prime[i] <= 0.0) {
	  ssa_prime[i] = 0.0;
	  g_prime[i] = 0.0;
	}
	else {
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
      for (i = 0; i < n; i++) {
	ext_prime[i] = ext[i];
	ssa_prime[i] = ssa[i];
	g_prime[i] = g[i];
      }
    }
    
    /* Calculate the source terms for the wide-angle calculation */
    for (i = 0; i < n; i++) {
      areal od = ext_prime[i]*drange; /* optical depth of layer */
      areal src_power;
      /* Fixed on 23 Feb 2012 */
      /*
      src_power = config->coherent_backscatter_enhancement
	*config->ss_multiplier*exp(-optical_depth)
	*(1.0-exp(-ext_prime[i]*ssa_prime[i]*drange));
	*/
      src_power = config->coherent_backscatter_enhancement
	*config->ss_multiplier*ssa_prime[i]*exp(-optical_depth)
	*(1.0-exp(-od));

      if (od > 0.0) {
	optical_depth += od;
      }

      if (g_prime[i] >= 1.0/(3.0*MS_MU1)) {
	/* If g exceeds 2/3, all energy is sent into outward diffuse stream */
	src_power_in[i] = 0.0;
	src_power_out[i] = src_power;
      }
      else if (g_prime[i] <= -1.0/(3.0*MS_MU1)) {
	/* If g is less than -2/3, all energy is sent into inward
	   diffuse stream (this never happens in atmospheric
	   scattering) */
	src_power_in[i] = src_power;
	src_power_out[i] = 0.0;
      }
      else {
	/* Two-stream phase function */
	src_power_in[i] = src_power*0.5*(1.0-3.0*g_prime[i]*MS_MU1);
	src_power_out[i] = src_power*0.5*(1.0+3.0*g_prime[i]*MS_MU1);
      }
      src_width2[i] = instrument.rho_transmitter*instrument.rho_transmitter
	*(range[i]-instrument.altitude)*(range[i]-instrument.altitude);
    }
  }

  else {
    areal transmittance = 1.0;
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
    
    for (i = 0; i < n; i++) {
      areal od, od_factor, power_factor, src_power, g_equiv, wide_factor;
      if (config->options & MS_SSA_SCALES_FORWARD_LOBE) {
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
	/* For single-scattering albedo of greater than 0.5/g, assume
	   all the absorption affects only the wide-angle part of the
	   phase function, not the forward lobe */
	/* If the extinction is zero then ssa, g, od_factor and
	   power_factor don't matter so we may as well do the same
	   thing here even if ssa=g=0 */
	od_factor = 1.0;
	power_factor = 2.0*ssa[i]-1.0;
      }
      else {
	/* For single-scattering albedo less than 0.5,
	   "diffraction-scaled" value of g becomes negative - this
	   doesn't make sense */
#ifdef CPPAD
	ms_real ssa_tmp = Value(CppAD::Var2Par(ssa[i]));
	ms_real g_tmp = Value(CppAD::Var2Par(g[i]));
	fprintf(stderr, "Error (CPPAD): cannot simulate multiple scattering including a forward lobe when ssa < 0.5/g and ext > 0 (ssa=%g, g=%g)\n", ssa_tmp, g_tmp);
#else
	fprintf(stderr, "Error: cannot simulate multiple scattering including a forward lobe when ssa < 0.5/g and ext > 0 (ssa=%g, g=%g)\n", value(ssa[i]), value(g[i]));
#endif

	return MS_FAILURE;
	//	od_factor = 2.0*(1.0-ssa[i]);
	//	power_factor = 0.0;
      }
      
      /* Calculate the optical depth to the start of gate i */
      if (ext_air) {
	od = (0.5*ext[i]*od_factor + ext_air[i])*drange;
      }
      else {
	od = 0.5*ext[i]*od_factor*drange;
      }
      /* Calculate the transmittance to gate i, including in-gate
	 attenuation */
      if (od > 0.0) {
	transmittance = exp(-optical_depth)*(1-exp(-od))/od;
	optical_depth += od;
	if (ext_air && ssa_air) {
	  wide_factor = (ext[i]*ssa[i] + ext_air[i]*ssa_air[i])
	    /(0.5*ext[i]*power_factor + ext_air[i]*ssa_air[i]);
	}
	else {
	  wide_factor = ext[i]*ssa[i]/(0.5*ext[i]*power_factor);
	}
      }
      else {
	transmittance = exp(-optical_depth);
	wide_factor = 1.0;
      }
      
      /* Calculate source power and source width-squared for the
	 wide-angle calculation */
      if (ext_air && ssa_air) {
	src_power = config->coherent_backscatter_enhancement
	  *config->ss_multiplier*transmittance
	  *(0.5*ext[i]*power_factor + ext_air[i]*ssa_air[i])*drange;
      }
      else {
	src_power = config->coherent_backscatter_enhancement
	  *config->ss_multiplier*transmittance
	  *0.5*ext[i]*power_factor*drange;
      }
      g_equiv = 1.0-wide_factor+wide_factor*g[i];
      if (g_equiv >= 1.0/(3.0*MS_MU1)) {
	src_power_in[i] = 0.0;
	src_power_out[i] = src_power;
      }
      else {
	src_power_in[i] = src_power*0.5*(1.0-3.0*g_equiv*MS_MU1);
	src_power_out[i] = src_power*0.5*(1.0+3.0*g_equiv*MS_MU1);
      }
    } /* End of spatial loop */
    
    if (ext_air && ssa_air) {
      /* Weighted delta scaling assuming that the gaseous asymmetry
	 factor is 0.0 */
      for (i = 0; i < n; i++) {
	ext_prime[i] = ext_air[i] + ext[i]*(1.0-ssa[i]*g[i]*g[i]);
	if (ext_air[i] + ext[i] <= 0.0) {
	  ssa_prime[i] = 0.0;
	  g_prime[i] = 0.0;
	}
	else {
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
      /* Standard delta scaling without air contribution */
      for (i = 0; i < n; i++) {
	ext_prime[i] = ext[i]*(1.0-ssa[i]*g[i]*g[i]);
	ssa_prime[i] = ssa[i]*(1.0-g[i]*g[i])/(1.0-ssa[i]*g[i]*g[i]);
	g_prime[i] = g[i]/(1.0+g[i]);
      }
    }
  }

  /*
  fprintf(stderr, "ext'  ssa'   g'   src_power_in   src_power_out   src_width2\n");
  for (i = 0; i < n; i++) {
    fprintf(stderr, "%10g %10g %10g %10g %10g %10g\n", ext_prime[i],
	    ssa_prime[i], g_prime[i], src_power_in[i], src_power_out[i],
	    src_width2[i]);
  }
  */


  /* Perform the 2-stream wide-angle multiple scattering
     calculation */
  //  if (first_gate == 0) {
    return ms_tdts(n, m, config, instrument, surface, 
		   range, ext_prime, ssa_prime, g_prime,
		   src_power_in, src_power_out, src_width2,
		   bscat_out);
    /*
  }
  else {
    return ms_tdts(n-first_gate, m-first_gate, config,
		   instrument, surface, 
		   range+first_gate, ext_prime+first_gate,
		   ssa_prime+first_gate, g_prime+first_gate,
		   src_power_in+first_gate, 
		   src_power_out+first_gate,
		   src_width2+first_gate,
		   bscat_out+first_gate);
  }
    */
}
