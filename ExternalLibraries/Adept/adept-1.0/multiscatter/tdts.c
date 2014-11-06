/* tdts.c -- Wide-angle multiple scattering, time-dependent two-stream method

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

#define ONE_OVER_FOUR_PI 0.0795774715459477

/* Simple calculation of transfer coefficients using an upwind Euler
   scheme. Note that this is unstable when the optical depth within a
   single layer approaches unity, i.e. when there is more than one
   scattering event in a single timestep. */
static
int
calculate_deltas_simple(
	    /* Inputs */
	    int n, const areal *ext, const areal *ssa, 
	    const areal *g, ms_real drange,
	    /* Outputs */
	    areal* delta0, areal* delta1, areal* delta2,
	    areal* delta3, areal* delta4, areal* delta5,
	    areal* transmittance, areal* transport_mfp)
{
  int i;
  areal optical_depth = 0.0;
  transmittance[-1] = 1.0;
  for (i = 0; i < n; i++) {
    areal gamma1 = (1-ssa[i]*0.5*(1+g[i]))/MS_MU1;
    areal gamma2 = ssa[i]*0.5*(1-g[i])/MS_MU1;
    delta0[i] = 1.0 - MS_MU1 - MS_MU1*drange*ext[i]*gamma1;
    delta1[i] = MS_MU1*drange*ext[i]*gamma2;
    delta2[i] = MS_MU1;
    delta3[i] = 0.0;
    delta4[i] = 0.0;
    delta5[i] = 0.0;
    if (ext[i] > 0.0) {
      transmittance[i] = exp(-optical_depth)
	*(1-exp(-ext[i]*drange))/(ext[i]*drange);
      optical_depth = optical_depth + ext[i]*drange;
      transport_mfp[i] = 1.0/((1.0-ssa[i]*g[i])*ext[i]);
    }
    else {
      transmittance[i] = exp(-optical_depth);
      transport_mfp[i] = 1.0e6*drange;
    }
  }
  delta1[-1] = delta1[n] = MS_MU1;
  delta3[-1] = delta3[n] = 0.0;
  delta4[-1] = delta5[n] = 0.0;
  transmittance[n] = transmittance[n-1];
  transport_mfp[-1] = transport_mfp[n] = 1.0e6*drange;
  return MS_SUCCESS;
}

/* Calculate transfer coefficients in a way that is stable even when
   the optical depth within a single layer approaches unity, i.e. when
   there is more than one scattering event in a single timestep. */
static
int
calculate_deltas(
	    /* Inputs */
	    int n, const areal *ext, const areal *ssa,
	    const areal *g, ms_real drange,
	    /* Outputs */
	    areal* delta0, areal* delta1, areal* delta2,
	    areal* delta3, areal* delta4, areal* delta5,
	    areal* transmittance, areal* transport_mfp)
{
  int i;
  ms_real dtrange = drange;
  areal optical_depth = 0.0;
  int negative_deltas = 0;
  int negative_deltas_after_mods = 0;
  transmittance[-1] = 1.0;
  for (i = 0; i < n; i++) {
    if (ext[i] > 0.0) {
      /* Mean free paths for transport and absorption */
      areal mfp_trans = 1.0/(ext[i]*(1-ssa[i]*g[i]));
      //      areal mfp_abs = 1.0/(ext[i]*(1-ssa[i]));
      /* Factors */
      areal f_trans = exp(-dtrange/mfp_trans);
      //      areal f_abs = exp(-dtrange/mfp_abs);
      areal f_abs = exp(-dtrange*ext[i]*(1.0-ssa[i]));
      areal f = MS_MU1*(mfp_trans/dtrange-f_trans/(1-f_trans));
      /* Diffusion terms */
      /* NEW: ssa[i] added 3 Aug 2011 */
      areal diffusion = MS_MU1*ssa[i]*f_abs*sqrt(mfp_trans/(3*dtrange));
      areal diffusion1 = diffusion*exp(-3.7*pow(mfp_trans/dtrange,0.75));
      areal diffusion2 = diffusion*exp(-3.7*mfp_trans/dtrange);
      /* Transfer terms */
      delta0[i] = f_trans*(1-MS_MU1) + (f_abs-f_trans)*(0.5-f)
	- diffusion1;
      delta1[i] = (f_abs-f_trans)*(1-f)*0.5
	- diffusion2;
      delta2[i] = MS_MU1*f_trans+(f_abs-f_trans)*f
	+ diffusion1*0.5;
      delta3[i] = (f_abs-f_trans)*f*0.25
	+ diffusion2*0.5;
      delta4[i] =
	+ diffusion1*0.5;
      delta5[i] = delta3[i];
      /* In principle, small negative terms are possible due to the
	 diffusion term exceeding the transport term: remove them but
	 conserving energy */
      if (delta0[i] < 0.0) {
	delta2[i] -= 0.5*delta0[i];
	delta4[i] -= 0.5*delta0[i];
	delta0[i] = 0.0;
	negative_deltas = 1;
      }
      if (delta1[i] < 0.0) {
	delta3[i] -= 0.5*delta1[i];
	delta5[i] -= 0.5*delta1[i];
	delta1[i] = 0.0;
	negative_deltas = 1;
      }
      if (delta2[i] < 0.0 || delta3[i] < 0.0 || delta4[i] < 0.0) {
	negative_deltas_after_mods = 1;
      }

      /* Transmittance and diffusivity */
      transmittance[i] = exp(-optical_depth)
	*(1-exp(-ext[i]*drange))/(ext[i]*drange);
      optical_depth = optical_depth + ext[i]*drange;      
      transport_mfp[i] = 1.0/((1.0-ssa[i]*g[i])*ext[i]);
    }
    else {
      delta0[i] = 1.0 - MS_MU1;
      delta1[i] = 0.0;
      delta2[i] = MS_MU1;
      delta3[i] = 0.0;
      delta4[i] = 0.0;
      delta5[i] = 0.0;
      transmittance[i] = exp(-optical_depth);
      transport_mfp[i] = 1.0e6*dtrange;
    }
  }
  delta1[-1] = delta1[n] = MS_MU1;
  delta3[-1] = delta3[n] = 0.0;
  delta4[-1] = delta5[n] = 0.0;
  transmittance[n] = transmittance[n-1];
  transport_mfp[-1] = transport_mfp[n] = 1.0e6*dtrange;
  if (negative_deltas) {
    fprintf(stderr, "Warning: delta0 or delta1 in the TDTS algorithm was negative: correction applied\n");
  }
  if (negative_deltas_after_mods) {
    fprintf(stderr, "Warning: at least one of delta2-delta5 in TDTS algorithm was negative: no correction applied\n");
  }
  return MS_SUCCESS;
}

/* The old method to estimate the lateral expansion of the photon
   distribution */

#ifdef COMPILE_OLD_EXPANSION_METHODS
static
areal
get_expansion_old(ms_real dt, areal D,
		  areal varw0, areal I, areal V)
{
  ms_real vmax = MS_C*sqrt(1.0-MS_MU1*MS_MU1);
#define W_FACTOR 1.0
#define DMAX_FACTOR 1.0
  areal Dmax, expansion;
  if (I <= 0.0) {
    Dmax = 1.0;
  }
  else if (V <= 0.0 || V <= I*varw0) {
    Dmax = DMAX_FACTOR*vmax*vmax*dt;
  }
  else {
    Dmax = DMAX_FACTOR*vmax*(vmax*dt + 2.0*sqrt(V/I - W_FACTOR*varw0));
  }
  expansion = dt*pow((1/(D*D*D) + 1/(Dmax*Dmax*Dmax)), -0.333333); 
  return expansion;
}


/* Estimate the lateral expansion of the photon distribution using a
   modified version of diffusion theory that accounts for the initial
   ballistic behaviour of photons. */

static
areal
get_expansion(ms_real drange, areal lt,
	      areal varw0, areal I, areal V)
{
  areal expansion, n;
  if (I <= 0.0) {
    return 0.0;
  }
  else if (V <= I*varw0) {
    n = 0.0;
  }
  else {
    areal norm_varw = (V/I-varw0)/(lt*lt);
    if (norm_varw < 0.8) {
      n = sqrt(norm_varw*2.0*MS_ONE_OVER_D_FACTOR);
    }
    else {
      n = 1.0 + norm_varw*MS_ONE_OVER_D_FACTOR;
    }
    if (norm_varw > 1.0e-3 && norm_varw < 4.0) {
      areal exp_minus_n = exp(-n);
      areal ln_predvarw_over_varw = log(MS_D_FACTOR*(n+exp_minus_n-1.0)
					  /norm_varw);
      areal ln_grad = (n+exp_minus_n-1.0)/(n-n*exp_minus_n);
      n *= exp(-ln_predvarw_over_varw*ln_grad);
    }
  }
  expansion = MS_D_FACTOR*lt*(drange + lt*(exp(-n-drange/lt)-exp(-n)));
  return expansion;
}
#endif

/* Estimate the lateral expansion of the photon distribution using a
   modified version of diffusion theory that accounts for the initial
   ballistic behaviour of photons: this is a faster version that is
   accurate to 1.5% */
static
//areal
void
get_expansion_fast(ms_real drange, areal lt,
		   areal varw0, areal I, areal V, areal* expansion)
{
  areal n;
  if (I <= 0.0) {
    //    return 0.0;  
    *expansion = 0.0;
    return;
  }
  else if (V <= I*varw0) {
    n = 0.0;
  }
  else {
    areal norm_varw = (V/I-varw0)/(lt*lt);
    if (norm_varw < 1.0) {
      areal sqrt_var = sqrt(norm_varw*2.0*MS_ONE_OVER_D_FACTOR);
      n = sqrt_var*(1.0+0.2*sqrt_var);
    }
    else {
      n = norm_varw*MS_ONE_OVER_D_FACTOR+1.0
	   -0.24/(norm_varw*norm_varw);
    }
  }
  *expansion = MS_D_FACTOR*lt*(drange + lt*(exp(-n-drange/lt)-exp(-n)));
  return;
}


/* Perform the Hogan and Battaglia (2008) time-dependent two-stream
   calculation for wide-angle multiple scattering, and add the result
   to bscat_out, to allow for another algorithm having previously
   calculated the single or QSA scattering return. */
int
ms_tdts(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const areal *ext,       /* Total extinction coefficient, m-1 */
    const areal *ssa,       /* Total single-scatter albedo */
    const areal *g,         /* Total asymmetry factor */
    const areal *src_power_in,/* Source function power: inwards */
    const areal *src_power_out,/* Source function power: outwards */
    const areal *src_width2,/* 1/e source function width^2, radians^2 */
    /* Output data */
    areal *bscat_out)       /* Measured backscatter, m-1 sr-1 */
{
  ms_real drange = fabs(range[2]-range[1]);
  
  int nt = 2*m;
  int it;
  int i;
  int ifov;
  
  /* Allocate the required vectors on the stack */
  areal _data[(n+2)*16];
  /* Incoming and outgoing radiances at the current and future
     timestep */
  areal* Iin              = _data+1;
  areal* Iout             = _data+(n+2)+1;
  areal* Iin_next         = _data+(n+2)*2+1;
  areal* Iout_next        = _data+(n+2)*3+1;
  /* Incoming and outgoing weighted variances at the current and future
     timestep */
  areal* Vin              = _data+(n+2)*4+1;
  areal* Vout             = _data+(n+2)*5+1;
  areal* Vin_next         = _data+(n+2)*6+1;
  areal* Vout_next        = _data+(n+2)*7+1;
  /* Transfer coefficients */
  areal* delta0           = _data+(n+2)*8+1;
  areal* delta1           = _data+(n+2)*9+1;
  areal* delta2           = _data+(n+2)*10+1;
  areal* delta3           = _data+(n+2)*11+1;
  areal* delta4           = _data+(n+2)*12+1;
  areal* delta5           = _data+(n+2)*13+1;
  /* Properties of the medium */
  areal* transport_mfp    = _data+(n+2)*14+1;
  areal* transmittance    = _data+(n+2)*15+1;

  /* Set the contents to zero */
  for (i = -1; i <= n; i++) {
    Iin[i] = 0.0;
    Iout[i] = 0.0;
    Iin_next[i] = 0.0;
    Iout_next[i] = 0.0;
    Vin[i] = 0.0;
    Vout[i] = 0.0;
    Vin_next[i] = 0.0;
    Vout_next[i] = 0.0;
    delta0[i] = 0.0;
    delta1[i] = 0.0;
    delta2[i] = 0.0;
    delta3[i] = 0.0;
    delta4[i] = 0.0;
    delta5[i] = 0.0;
    transport_mfp[i] = 0.0;
    transmittance[i] = 0.0;
  }

  /* Calculate transfer coefficients */
  if (config->options & MS_SIMPLE_2S_COEFFTS) {
    calculate_deltas_simple(n, ext, ssa, g, drange,
			    delta0, delta1, delta2,
			    delta3, delta4, delta5,
			    transmittance, transport_mfp);
  }
  else {
    calculate_deltas(n, ext, ssa, g, drange,
		     delta0, delta1, delta2,
		     delta3, delta4, delta5,
		     transmittance, transport_mfp);
  }

  config->total_src = 0.0;
  config->total_reflected = 0.0;

  /* Loop through each timestep */
  for (it = 0; it < nt; it++) {
    areal *tmp_I;
    areal expansion;
    /* Work out how many spatial points to calculate given that we
       don't want to waste time simulating regions not yet reached by
       the outgoing beam, or regions that can't scatter back to the
       receiver within the time period of interest */
    int max_i = it;
    if (max_i >= n) {
      max_i = n-1;
    }
    /* If the internal radiances and variances are to be output,
       then we calculate even in regions that will not scatter
       back to the receiver within the sample time. */
    if (!(config->options & MS_PROPAGATION_TO_STDERR) && (it > nt-n)) {
      max_i = nt-it;
    }

    config->total_reflected += value(delta1[0]*Iin[0]);

    /* Loop through each range gate */
    for (i = 0; i <= max_i ; i++) {

      /* Step the photon energy forward in time */
      Iin_next[i]
	= Iin[i]   * delta0[i]
	+ Iout[i]  * delta1[i]
	+ Iin[i+1] * delta2[i+1]
	+ Iout[i-1]* delta3[i-1]
	+ Iin[i-1] * delta4[i-1]
	+ Iout[i+1]* delta5[i+1];
      Iout_next[i]
	= Iout[i]  * delta0[i]
	+ Iin[i]   * delta1[i]
	+ Iout[i-1]* delta2[i-1]
	+ Iin[i+1] * delta3[i+1]
	+ Iout[i+1]* delta4[i+1]
	+ Iin[i-1] * delta5[i-1];
      if (Iin_next[i] < 0.0) {
	fprintf(stderr, "Warning: negative radiance at line %d of %s, timestep %d gate %d\n"
		"   (Iin_next=%g Iin=%g Iout=%g delta0=%g delta1=%g delta2=%g)\n",
		__LINE__, __FILE__, it, i,
		value(Iin_next[i]), value(Iin[i]), value(Iout[i]),
		value(delta0[i]), value(delta1[i+1]), value(delta2[i]));
      }

      /* Step the variance forward in time */
      /* #define USE_SLOW_EXPANSION 1 */
#ifdef USE_SLOW_EXPANSION
      expansion = get_expansion(drange, transport_mfp[i], src_width2[i],
				Iin[i], Vin[i]);
#else
      get_expansion_fast(drange, transport_mfp[i],
			 src_width2[i],
			 Iin[i], Vin[i], &expansion);
#endif
      Vin_next[i]
	= Vin[i]   * delta0[i]
	+ Vout[i]  * delta1[i]
	+ Vin[i+1] * delta2[i+1]
	+ Vout[i-1]* delta3[i-1]
	+ Vin[i-1] * delta4[i-1]
	+ Vout[i+1]* delta5[i+1]
	+ Iin[i]   * expansion;
#ifdef USE_SLOW_EXPANSION
      expansion = get_expansion(drange, transport_mfp[i], src_width2[i],
				Iout[i], Vout[i]);
#else
      get_expansion_fast(drange, transport_mfp[i],
			 src_width2[i],
			 Iout[i], Vout[i], &expansion);
#endif
      Vout_next[i]
	= Vout[i]  * delta0[i]
	+ Vin[i]   * delta1[i]
	+ Vout[i-1]* delta2[i-1]
	+ Vin[i+1] * delta3[i+1]
	+ Vout[i+1]* delta4[i+1]
	+ Vin[i-1] * delta5[i-1]
	+ Iout[i]  * expansion;

      /* Increment the apparent backscatter, looping over each
	 field-of-view */
      for (ifov = 0; ifov < instrument.nfov; ifov++) {
	ms_real footprint_radius2 
	  = instrument.rho_receiver[ifov]*instrument.rho_receiver[ifov]
	  *(range[i]-instrument.altitude)*(range[i]-instrument.altitude);
	int ibscat = ifov*m + (it+i)/2; /* BEST FOR I3RC */
	/* Check that we are not going to overwrite the next point in
	   memory */
	if (ibscat >= (ifov+1)*m) continue;
	/*      ibscat = ifov*n + (it+i-1)/2; // BEST FOR RADAR SCENARIO 3 */
	if (Vin_next[i] > 0.0 && Vout_next[i] > 0.0) {
	  areal bscat_inc;
	  if (instrument.receiver_type == MS_TOP_HAT) {
	    /* Lidar-type top-hat receiver */
	    /* First restrict superluminal travel */
	    ms_real footprint_radius_max = (it-i) * drange;
	    ms_real footprint_radius_max2 
	      = footprint_radius_max*footprint_radius_max 
	      + instrument.rho_transmitter*instrument.rho_transmitter
	      *(range[i]-instrument.altitude)*(range[i]-instrument.altitude);
	    ms_real eff_footprint_radius2 = footprint_radius2;
	    if (footprint_radius_max2 < footprint_radius2) {
	      eff_footprint_radius2 = footprint_radius_max2;
	    }
	    if (g[i] >= 1.0/(3.0*MS_MU1)) {
	      bscat_inc = transmittance[i]
		*ONE_OVER_FOUR_PI
		*(Iin_next[i] * ssa[i]*ext[i]*2.0
		  *(1.0-exp(-(eff_footprint_radius2*Iin_next[i]/Vin_next[i]))));
	    }
	    else {
	      bscat_inc = transmittance[i]
		*ONE_OVER_FOUR_PI
		*(Iin_next[i] * ssa[i]*ext[i]*(1.0+3.0*g[i]*MS_MU1)
		  *(1.0-exp(-(eff_footprint_radius2*Iin_next[i]/Vin_next[i])))
		  +Iout_next[i] * ssa[i]*ext[i]*(1.0-3.0*g[i]*MS_MU1)
		  *(1.0-exp(-(eff_footprint_radius2*Iout_next[i]/Vout_next[i]))));
	    }
	  }
	  else {
	    /* Radar-type Gaussian receiver */
	    if (g[i] >= 1.0/(3.0*MS_MU1)) {
	      bscat_inc = transmittance[i]
		*ONE_OVER_FOUR_PI
		*(Iin_next[i] * ssa[i]*ext[i]*2.0
		  /(1.0+Vin_next[i]/(Iin_next[i]*footprint_radius2)));
	    }
	    else {
	      bscat_inc = transmittance[i]
		*ONE_OVER_FOUR_PI
		*(Iin_next[i] * ssa[i]*ext[i]*(1.0+3.0*g[i]*MS_MU1)
		  /(1.0+Vin_next[i]/(Iin_next[i]*footprint_radius2))
		  +Iout_next[i] * ssa[i]*ext[i]*(1.0-3.0*g[i]*MS_MU1)
		  /(1.0+Vout_next[i]/(Iout_next[i]*footprint_radius2)));
	    }
	  }

	  /* Assign backscatter increment */
	  bscat_out[ibscat] += bscat_inc;
	  if (config->options & MS_ANNULAR_DETECTORS) {
	    if (ifov < instrument.nfov-1) {
	      /* This backscatter should be subtracted from the next
		 field-of-view */
	      bscat_out[ibscat+m] -= bscat_inc;
	    }
	  }
	}
      }
    }
    
    /* Print propagation variables (radiances and variances) to
       standard error */
    if (config->options & MS_PROPAGATION_TO_STDERR) {
      for (i = 1; i < n ; i++) {
	/*      for (i = 0; i < n ; i++) { */
	if (it < n && i == it) {
	  fprintf(stderr, "%d %d %g %g %g %g %g %g %g %g\n", it, i,
		  value(Iout_next[i]), value(Iin_next[i]),
		  value(Vout_next[i]), value(Vin_next[i]),
		  value(src_power_out[it]), value(src_power_in[it]),
		  value(src_power_out[it])*value(src_width2[it]),
		  value(src_power_in[it])*value(src_width2[it]));
	}
	else {
	  fprintf(stderr, "%d %d %g %g %g %g 0 0 0 0\n", it, i,
		  value(Iout_next[i]), value(Iin_next[i]),
		  value(Vout_next[i]), value(Vin_next[i]));
	}
      }
    }

    /* Add the incoming source power if we are in the first half of
       the time simulated */
    if (it < n) {
      /* The source power is split into that coming in towards the
	 instrument and that travelling out */
      if (src_power_in[it] < 0.0 || src_power_out[it] < 0.0) {
	fprintf(stderr, "Error: negative source terms for two-stream equations at line %d of %s\n"
		        "       (src_power_out[%d]=%g src_power_in[%d]=%g ssa=%g g=%g)\n",
		__LINE__, __FILE__,
		it, value(src_power_out[it]), it, value(src_power_in[it]), value(ssa[it]), value(g[it]));
	return MS_FAILURE;
      }
      else if (src_width2[it] < 0.0) {
	fprintf(stderr, "Error: negative source width for two-stream equations at line %d of %s\n",
		__LINE__, __FILE__);
	return MS_FAILURE;
      }
      Iin_next[it] += src_power_in[it];
      Iout_next[it] += src_power_out[it];
      Vin_next[it] += src_power_in[it] * src_width2[it];
      Vout_next[it] += src_power_out[it] * src_width2[it];
      config->total_src += value(src_power_in[it] + src_power_out[it]);
    }

#define swap_I(I, I_next) tmp_I=I; I=I_next; I_next=tmp_I;
    swap_I(Iin, Iin_next);
    swap_I(Iout, Iout_next);
    swap_I(Vin, Vin_next);
    swap_I(Vout, Vout_next);
  }

  return MS_SUCCESS;
}
