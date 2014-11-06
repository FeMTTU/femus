/* tdts_AD.c -- Adjoint of wide-angle multiple scattering, time-dependent two-stream method

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

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ms.h"

#define ONE_OVER_FOUR_PI 0.0795774715459477
#define MINIMUM_OD 1.0e-8

/* Calculate transfer coefficients in a way that is stable even when
   the optical depth within a single layer approaches unity, i.e. when
   there is more than one scattering event in a single timestep. */
static
int
calculate_deltas(
	    /* Inputs */
	    int n, const ms_real *ext, const ms_real *ssa,
	    const ms_real *g, ms_real drange,
	    /* Outputs */
	    ms_real* delta0, ms_real* delta1, ms_real* delta2,
	    ms_real* delta3, ms_real* delta4, ms_real* delta5,
	    ms_real* transmittance, ms_real* transport_mfp)
{
  int i;
  ms_real dtrange = drange;
  ms_real optical_depth = 0.0;
  transmittance[-1] = 1.0;
  for (i = 0; i < n; i++) {
    if (ext[i] > 0.0) {
      /* Mean free paths for transport and absorption */
      ms_real mfp_trans = 1.0/(ext[i]*(1.0-ssa[i]*g[i]));
      ms_real mfp_abs = 1.0/(ext[i]*(1.0-ssa[i]));
      /* Factors */
      ms_real f_trans = exp(-dtrange/mfp_trans);
      ms_real f_abs = exp(-dtrange/mfp_abs);
      ms_real f = MS_MU1*(mfp_trans/dtrange-f_trans/(1.0-f_trans));
      /* Diffusion terms */
      /* NEW: ssa[i] added 3 Aug 2011 */
      ms_real diffusion = MS_MU1*ssa[i]*f_abs*sqrt(mfp_trans/(3*dtrange));
      ms_real diffusion1 = diffusion*exp(-3.7*pow(mfp_trans/dtrange,0.75));
      ms_real diffusion2 = diffusion*exp(-3.7*mfp_trans/dtrange);
      /* Transfer terms */
      delta0[i] = f_trans*(1.0-MS_MU1) + (f_abs-f_trans)*(0.5-f)
	- diffusion1;
      delta1[i] = (f_abs-f_trans)*(1.0-f)*0.5
	- diffusion2;
      delta2[i] = MS_MU1*f_trans+(f_abs-f_trans)*f
	+ diffusion1*0.5;
      delta3[i] = (f_abs-f_trans)*f*0.25
	+ diffusion2*0.5;
      delta4[i] =
	+ diffusion1*0.5;
      delta5[i] = delta3[i];
      if (0 && i == 0) {
	delta2[i] += delta5[i];
	delta5[i] = 0.0;
	delta1[i] /= 2.0;
      }
      /* Transmittance and diffusivity */
      transmittance[i] = exp(-optical_depth)
	*(1.0-exp(-ext[i]*drange))/(ext[i]*drange);
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
      transport_mfp[i] = dtrange/MINIMUM_OD;
    }
  }
  delta1[-1] = delta1[n] = MS_MU1;
  delta3[-1] = delta3[n] = 0.0;
  delta4[-1] = delta5[n] = 0.0;
  transmittance[n] = transmittance[n-1];
  transport_mfp[-1] = transport_mfp[n] = dtrange/MINIMUM_OD;
  return MS_SUCCESS;
}

/* Calculate transfer coefficients in a way that is stable even when
   the optical depth within a single layer approaches unity, i.e. when
   there is more than one scattering event in a single timestep. */
static
int
calculate_deltas_AD(
	    /* Inputs */
	    int n, 
	    const ms_real* ext, const ms_real* ssa, const ms_real* g, 
	    ms_real drange,
	    /* Outputs */
	    ms_real* delta0, ms_real* delta1, ms_real* delta2,
	    ms_real* delta3, ms_real* delta4, ms_real* delta5,
	    ms_real* transmittance, ms_real* transport_mfp,
	    /* Adjoint inputs/outputs */
	    ms_real* delta0_AD, ms_real* delta1_AD, 
	    ms_real* delta2_AD, ms_real* delta3_AD,
	    ms_real* delta4_AD, ms_real* delta5_AD,
	    ms_real* transmittance_AD, ms_real* transport_mfp_AD,
	    /* Adjoint outputs */
	    ms_real* ext_AD, ms_real* ssa_AD, ms_real* g_AD)
{
  int i;
  ms_real dtrange = drange;
  ms_real transmittance_mid[n+1]; /* transmittance to half-level */
  ms_real transmittance_mid_AD = 0.0;
  transmittance_mid[0] = transmittance[-1] = 1.0;
  for (i = 0; i < n; i++) {
    if (ext[i] > 0.0) {
      /* Mean free paths for transport and absorption */
      ms_real mfp_trans = 1.0/(ext[i]*(1.0-ssa[i]*g[i]));
      ms_real mfp_abs = 1.0/(ext[i]*(1.0-ssa[i]));
      /* Factors */
      ms_real f_trans = exp(-dtrange/mfp_trans);
      ms_real f_abs = exp(-dtrange/mfp_abs);
      ms_real f = MS_MU1*(mfp_trans/dtrange-f_trans/(1.0-f_trans));
      /* Diffusion terms */
      /* NEW: ssa[i] added 3 Aug 2011 */
      ms_real diffusion = MS_MU1*ssa[i]*f_abs*sqrt(mfp_trans/(3*dtrange));
      ms_real diffusion1 = diffusion*exp(-3.7*pow(mfp_trans/dtrange,0.75));
      ms_real diffusion2 = diffusion*exp(-3.7*mfp_trans/dtrange);
      /* Transfer terms */
      delta0[i] = f_trans*(1.0-MS_MU1) + (f_abs-f_trans)*(0.5-f)
	- diffusion1;
      delta1[i] = (f_abs-f_trans)*(1.0-f)*0.5
	- diffusion2;
      delta2[i] = MS_MU1*f_trans+(f_abs-f_trans)*f
	+ diffusion1*0.5;
      delta3[i] = (f_abs-f_trans)*f*0.25
	+ diffusion2*0.5;
      delta4[i] =
	+ diffusion1*0.5;
      delta5[i] = delta3[i];
      /* Transmittance and diffusivity */
      transmittance_mid[i+1]
	= transmittance_mid[i]*exp(-ext[i]*drange);
      transmittance[i] = (transmittance_mid[i]-transmittance_mid[i+1])
	/ (ext[i]*drange);

      transport_mfp[i] = 1.0/((1.0-ssa[i]*g[i])*ext[i]);
    }
    else {
      delta0[i] = 1.0 - MS_MU1;
      delta1[i] = 0.0;
      delta2[i] = MS_MU1;
      delta3[i] = 0.0;
      delta4[i] = 0.0;
      delta5[i] = 0.0;
      transmittance[i] = transmittance_mid[i];
      transmittance_mid[i+1] = transmittance_mid[i];
      transport_mfp[i] = dtrange/MINIMUM_OD;
    }
  }
  delta1[-1] = delta1[n] = MS_MU1;
  delta3[-1] = delta3[n] = 0.0;
  delta4[-1] = delta5[n] = 0.0;
  transmittance[n] = transmittance[n-1];
  transport_mfp[-1] = transport_mfp[n] = dtrange/MINIMUM_OD;

  /* Adjoint */
  delta1_AD[-1] = delta1_AD[n] = 0.0;
  delta3_AD[-1] = delta3_AD[n] = 0.0;
  delta4_AD[-1] = delta5_AD[n] = 0.0;
  transmittance_AD[n-1] += transmittance_AD[n];
  transmittance_AD[n] = 0.0;
  transport_mfp_AD[-1] = transport_mfp_AD[n] = 0.0;

  for (i = n-1; i >= 0; i--) {
    ms_real my_ext = ext[i];
    if (my_ext <= 0.0) {
      /* Minimum optical depth of layer is 1.0e-8, in order that
	 adjoint calculation works for ext[i] = 0 */
      my_ext = MINIMUM_OD/drange;
    }
    if (ext[i] > 0.0) {
      /* Define required terms */
      /* Mean free paths for transport and absorption */
      ms_real mfp_trans = 1.0/(my_ext*(1.0-ssa[i]*g[i]));
      /* Factors */
      ms_real f_trans = exp(-dtrange/mfp_trans);
      ms_real f_abs = exp(-dtrange*my_ext*(1.0-ssa[i]));
      ms_real f = MS_MU1*(mfp_trans/dtrange-f_trans/(1.0-f_trans));
      /* Diffusion terms */
      ms_real diffusion_sqrt = sqrt(mfp_trans/(3*dtrange));
      ms_real diffusion = MS_MU1*ssa[i]*f_abs*diffusion_sqrt;
      ms_real diffusion1_exp = exp(-3.7*pow(mfp_trans/dtrange,0.75));
      ms_real diffusion1 = diffusion*diffusion1_exp;
      ms_real diffusion2_exp = exp(-3.7*mfp_trans/dtrange);
      ms_real diffusion2 = diffusion*diffusion2_exp;
      /* Now do the adjoints */
      ms_real f_trans_AD = 0.0, f_abs_AD = 0.0, f_AD = 0.0, mfp_trans_AD = 0.0;
      ms_real diffusion1_AD = 0.0, diffusion2_AD = 0.0, diffusion_AD = 0.0;
      /* Transmittance and diffusivity */
      ms_real factor 
	= transport_mfp_AD[i]*transport_mfp[i] / (1.0-ssa[i]*g[i]);
      ssa_AD[i] += g[i]*factor;
      g_AD[i] += ssa[i]*factor;

      ext_AD[i] -= transport_mfp_AD[i]*transport_mfp[i]/my_ext;
      transport_mfp_AD[i] = 0.0;

      factor = transmittance_AD[i] / (my_ext*drange);
      ext_AD[i] -= (transmittance[i]*transmittance_AD[i]/my_ext
		    +drange*transmittance_mid[i+1]
		    *(transmittance_mid_AD-factor));
      transmittance_mid_AD = factor + exp(-my_ext*drange)
	*(transmittance_mid_AD - factor);

      transmittance_AD[i] = 0.0;

      /* Transfer terms */
      f_trans_AD = (0.5+f-MS_MU1)*delta0_AD[i] 
	+ 0.5*(f-1.0)*delta1_AD[i]
	+ (MS_MU1 - f)*delta2_AD[i]
	-0.25*f*(delta3_AD[i] + delta5_AD[i]);

      f_abs_AD = (0.5-f)*delta0_AD[i] 
	+ 0.5*(1.0-f)*delta1_AD[i]
	+ f*delta2_AD[i]
	+ 0.25*f*(delta3_AD[i] + delta5_AD[i]);

      f_AD = (f_trans-f_abs)*delta0_AD[i] 
	+ 0.5*(f_trans-f_abs)*delta1_AD[i]
	+ (f_abs-f_trans)*delta2_AD[i]
	+ 0.25*(f_abs-f_trans)*(delta3_AD[i]+delta5_AD[i]);

      diffusion1_AD = -delta0_AD[i] 
	+ 0.5*(delta2_AD[i] + delta4_AD[i]);

      diffusion2_AD = -delta1_AD[i]
	+ 0.5*(delta3_AD[i] + delta5_AD[i]);

      delta0_AD[i] = delta1_AD[i] = delta2_AD[i]
	= delta3_AD[i] = delta4_AD[i] = delta5_AD[i] = 0.0;

      /* Diffusion terms */
      diffusion_AD = diffusion2_exp*diffusion2_AD
	+ diffusion1_exp*diffusion1_AD;
      mfp_trans_AD = -3.7*diffusion2*diffusion2_AD/dtrange
	- 3.7*diffusion1*0.75*pow(mfp_trans,-0.25)
	   *pow(dtrange,-0.75)*diffusion1_AD;
      diffusion1_AD = diffusion2_AD = 0.0;

      ssa_AD[i] += MS_MU1*diffusion_sqrt*f_abs*diffusion_AD;
      f_abs_AD += MS_MU1*diffusion_sqrt*ssa[i]*diffusion_AD;
      mfp_trans_AD += MS_MU1*f_abs*ssa[i]*0.5*diffusion_sqrt*diffusion_AD/mfp_trans;
      diffusion_AD = 0.0;

      /* Factors */
      mfp_trans_AD += MS_MU1*f_AD / dtrange;
      f_trans_AD -= MS_MU1*(1.0+f_trans/(1.0-f_trans))*f_AD/(1.0-f_trans);
      f_AD = 0.0;

      ext_AD[i] -= f_abs*dtrange*(1.0-ssa[i])*f_abs_AD;
      ssa_AD[i] += f_abs*dtrange*my_ext*f_abs_AD;
      f_abs_AD = 0.0;

      mfp_trans_AD += f_trans*dtrange*f_trans_AD/(mfp_trans*mfp_trans);
      f_trans_AD = 0.0;

     
      /* Mean free paths for transport and absorption */
      ext_AD[i] -= mfp_trans*mfp_trans_AD/my_ext;
      ssa_AD[i] += mfp_trans*mfp_trans_AD*g[i]/(1.0-ssa[i]*g[i]);
      g_AD[i] += mfp_trans*mfp_trans_AD*ssa[i]/(1.0-ssa[i]*g[i]);
      mfp_trans_AD = 0.0;
      
    }
    else {
      delta0_AD[i] = delta1_AD[i] = delta2_AD[i]
	= delta3_AD[i] = delta4_AD[i] = delta5_AD[i] = 0.0;
      transmittance_mid_AD += transmittance_AD[i];
      transport_mfp_AD[i] = 0.0;
    }
      
  }

  transmittance_AD[-1] = 0.0;

  return MS_SUCCESS;
}


/* Estimate the lateral expansion of the photon distribution using a
   modified version of diffusion theory that accounts for the initial
   ballistic behaviour of photons: this is a faster version that is
   accurate to 1.5% */
static
ms_real
get_expansion_fast(ms_real drange, ms_real lt,
		   ms_real varw0, ms_real I, ms_real V)
{
  ms_real expansion, n;
  if (I <= 0.0) {
    return 0.0;
  }
  else if (V <= I*varw0) {
    n = 0.0;
  }
  else {
    ms_real norm_varw = (V/I-varw0)/(lt*lt);
    if (norm_varw < 1.0) {
      ms_real sqrt_var = sqrt(norm_varw*2.0*MS_ONE_OVER_D_FACTOR);
      n = sqrt_var*(1.0+0.2*sqrt_var);
    }
    else {
      n = norm_varw*MS_ONE_OVER_D_FACTOR+1.0
	   -0.24/(norm_varw*norm_varw);
    }
  }
  expansion = MS_D_FACTOR*lt*(drange + lt*(exp(-n-drange/lt)-exp(-n)));
  return expansion;
}

/* Adjoint version */
static
ms_real
get_expansion_fast_AD(ms_real drange, /* Parameters */
		      /* Inputs */
		      ms_real lt, ms_real varw0, 
		      ms_real I, ms_real V,
		      ms_real* expansion_AD, /* Adjoint input */
		      /* Adjoint outputs: */
		      ms_real* lt_AD, ms_real* varw0_AD,
		      ms_real* I_AD, ms_real* V_AD)
{
  ms_real expansion, n;
  ms_real norm_varw = 0.0;
  ms_real sqrt_var = 0.0;
  ms_real exp1, exp2;
  ms_real n_AD = 0.0, norm_varw_AD = 0.0;

  if (I <= 0.0) {
    return 0.0;
  }
  else if (V <= I*varw0) {
    n = 0.0;
  }
  else {
    norm_varw = (V/I-varw0)/(lt*lt);
    if (norm_varw < 1.0) {
      sqrt_var = sqrt(norm_varw*2.0*MS_ONE_OVER_D_FACTOR);
      n = sqrt_var*(1.0+0.2*sqrt_var);
    }
    else {
      n = norm_varw*MS_ONE_OVER_D_FACTOR+1.0
	   -0.24/(norm_varw*norm_varw);
    }
  }
  exp1 = exp(-n-drange/lt);
  exp2 = exp(-n);
  expansion = MS_D_FACTOR*lt*(drange + lt*(exp1-exp2));
  
  /* Adjoint part */
  *lt_AD += MS_D_FACTOR*(drange + 2*lt*(exp1-exp2)
			 +drange*exp1)* *expansion_AD;
  n_AD = MS_D_FACTOR*lt*lt*(exp2-exp1)* *expansion_AD;
  *expansion_AD = 0.0;

  if (V <= I*varw0) {
    n_AD = 0.0;
    return expansion;
  }
  else if (norm_varw < 1.0) {
    ms_real sqrt_var_AD = (1.0 + 0.4*sqrt_var)*n_AD;
    n_AD = 0.0;
    norm_varw_AD = sqrt(MS_ONE_OVER_D_FACTOR*0.5/norm_varw)*sqrt_var_AD;
    sqrt_var_AD = 0.0;
  }
  else {
    norm_varw_AD = n_AD*(MS_ONE_OVER_D_FACTOR 
			 + (2.0*0.24) / (norm_varw*norm_varw*norm_varw));
    n_AD = 0.0;
  }

  *V_AD += norm_varw_AD / (I*lt*lt);
  *I_AD -= norm_varw_AD * V / (I*I*lt*lt);
  if (varw0_AD) {
    /* Sometimes the adjoint of the variance not required */
    *varw0_AD -= norm_varw_AD / (lt*lt);
  }
  *lt_AD -= norm_varw_AD * 2.0 * norm_varw / lt;
  norm_varw_AD = 0.0;    

  return expansion;
}


/* Perform the Hogan and Battaglia (2008) time-dependent two-stream
   calculation for wide-angle multiple scattering, and add the result
   to bscat_out, to allow for another algorithm having previously
   calculated the single or QSA scattering return. */
int
ms_tdts_AD(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *ext,       /* Total extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *src_power_in,/* Source function power: inwards */
    const ms_real *src_power_out,/* Source function power: outwards */
    const ms_real *src_width2,/* 1/e source function width^2, radians^2 */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    /* Adjoint inputs */
    const ms_real *bscat_AD,
    /* Adjoint outputs */
    ms_real *ext_AD,
    ms_real *ssa_AD,
    ms_real *g_AD,
    ms_real *src_power_in_AD,
    ms_real *src_power_out_AD,
    ms_real *src_width2_AD)
{
  ms_real drange = fabs(range[2]-range[1]);
  
  int nt = 2*m;
  int it;
  int i;
  int ifov;
  
  /* Incoming and outgoing radiances at all timesteps */
  ms_real Iin[nt+1][n+2];
  ms_real Iout[nt+1][n+2];
  /* Incoming and outgoing weighted variances at all timesteps */
  ms_real Vin[nt+1][n+2];
  ms_real Vout[nt+1][n+2];
  /* Allocate the required vectors on the stack */
  ms_real _data[(n+2)*16];
  /* Transfer coefficients */
  ms_real* delta0           = _data+1;
  ms_real* delta1           = _data+(n+2)+1;
  ms_real* delta2           = _data+(n+2)*2+1;
  ms_real* delta3           = _data+(n+2)*3+1;
  ms_real* delta4           = _data+(n+2)*4+1;
  ms_real* delta5           = _data+(n+2)*5+1;
  /* Properties of the medium */
  ms_real* transport_mfp    = _data+(n+2)*6+1;
  ms_real* transmittance    = _data+(n+2)*7+1;
  /* Adjoints */
  ms_real Iin_AD[n+2];
  ms_real Iout_AD[n+2];
  ms_real Vin_AD[n+2];
  ms_real Vout_AD[n+2];
  ms_real Iin_next_AD[n+2];
  ms_real Iout_next_AD[n+2];
  ms_real Vin_next_AD[n+2];
  ms_real Vout_next_AD[n+2];
  ms_real* delta0_AD        = _data+(n+2)*8+1;
  ms_real* delta1_AD        = _data+(n+2)*9+1;
  ms_real* delta2_AD        = _data+(n+2)*10+1;
  ms_real* delta3_AD        = _data+(n+2)*11+1;
  ms_real* delta4_AD        = _data+(n+2)*12+1;
  ms_real* delta5_AD        = _data+(n+2)*13+1;
  ms_real* transport_mfp_AD = _data+(n+2)*14+1;
  ms_real* transmittance_AD = _data+(n+2)*15+1;


  /* Set the contents to zero */
  for (i = -1; i <= n; i++) {
    int i1 = i+1;
    int j;
    for (j = 0; j <= nt; j++) {
      Iin[j][i1] = 0.0;
      Iout[j][i1] = 0.0;
      Vin[j][i1] = 0.0;
      Vout[j][i1] = 0.0;
    }
    delta0[i] = 0.0;
    delta1[i] = 0.0;
    delta2[i] = 0.0;
    delta3[i] = 0.0;
    delta4[i] = 0.0;
    delta5[i] = 0.0;
    transport_mfp[i] = 0.0;
    transmittance[i] = 0.0;
    /* Adjoints */
    Iin_AD[i1] = 0.0;
    Iout_AD[i1] = 0.0;
    Vin_AD[i1] = 0.0;
    Vout_AD[i1] = 0.0;
    Iin_next_AD[i1] = 0.0;
    Iout_next_AD[i1] = 0.0;
    Vin_next_AD[i1] = 0.0;
    Vout_next_AD[i1] = 0.0;
    delta0_AD[i] = 0.0;
    delta1_AD[i] = 0.0;
    delta2_AD[i] = 0.0;
    delta3_AD[i] = 0.0;
    delta4_AD[i] = 0.0;
    delta5_AD[i] = 0.0;
    transport_mfp_AD[i] = 0.0;
    transmittance_AD[i] = 0.0;
  }

  /* Calculate transfer coefficients */
  calculate_deltas(n, ext, ssa, g, drange,
		   delta0, delta1, delta2,
		   delta3, delta4, delta5,
		   transmittance, transport_mfp);

  config->total_src = 0.0;
  config->total_reflected = 0.0;

  /* Loop through each timestep */
  for (it = 0; it < nt; it++) {
    ms_real expansion;
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

    /*    max_i = n-1; */

    Iin[it][0] = 0.0;
    Iout[it][0] = 0.0;
    Vin[it][0] = 0.0;
    Vout[it][0] = 0.0;

    /* Loop through each range gate */
    for (i = 0; i <= max_i ; i++) {
      /* Spatial index for Iin, Iout, Vin, Vout is offset by one: */
      int i1 = i+1;

      /* Step the photon energy forward in time */
      Iin[it+1][i1]
	= Iin[it][i1]   * delta0[i]
	+ Iout[it][i1]  * delta1[i]
	+ Iin[it][i1+1] * delta2[i+1]
	+ Iout[it][i1-1]* delta3[i-1]
	+ Iin[it][i1-1] * delta4[i-1]
	+ Iout[it][i1+1]* delta5[i+1];
      Iout[it+1][i1]
	= Iout[it][i1]  * delta0[i]
	+ Iin[it][i1]   * delta1[i]
	+ Iout[it][i1-1]* delta2[i-1]
	+ Iin[it][i1+1] * delta3[i+1]
	+ Iout[it][i1+1]* delta4[i+1]
	+ Iin[it][i1-1] * delta5[i-1];


      if (Iin[it+1][i1] < 0.0) {
	fprintf(stderr, "Warning: negative radiance at line %d of %s, timestep %d gate %d\n"
		"   (Iin_next=%g Iin=%g Iout=%g delta0=%g delta1=%g delta2=%g)\n",
		__LINE__, __FILE__, it, i,
		Iin[it+1][i1], Iin[it][i1], Iout[it][i1],
		delta0[i], delta1[i+1], delta2[i]);
      }

      /* Step the variance forward in time */
      expansion = get_expansion_fast(drange, transport_mfp[i],
				     src_width2[i],
				     Iin[it][i1], Vin[it][i1]);
      Vin[it+1][i1]
	= Vin[it][i1]   * delta0[i]
	+ Vout[it][i1]  * delta1[i]
	+ Vin[it][i1+1] * delta2[i+1]
	+ Vout[it][i1-1]* delta3[i-1]
	+ Vin[it][i1-1] * delta4[i-1]
	+ Vout[it][i1+1]* delta5[i+1]
	+ Iin[it][i1]   * expansion;


      expansion = get_expansion_fast(drange, transport_mfp[i],
				     src_width2[i],
				     Iout[it][i1], Vout[it][i1]);
      Vout[it+1][i1]
	= Vout[it][i1]  * delta0[i]
	+ Vin[it][i1]   * delta1[i]
	+ Vout[it][i1-1]* delta2[i-1]
	+ Vin[it][i1+1] * delta3[i+1]
	+ Vout[it][i1+1]* delta4[i+1]
	+ Vin[it][i1-1] * delta5[i-1]
	+ Iout[it][i1]  * expansion;

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
	if (Vin[it+1][i1] > 0.0 && Vout[it+1][i1] > 0.0) {
	  ms_real bscat_inc;
	  if (instrument.receiver_type == MS_TOP_HAT) {
	    /* Lidar-type top-hat receiver */
	    /* First restrict superluminal travel */
	    ms_real footprint_radius_max = (it-i) * drange;
	    ms_real footprint_radius_max2 = footprint_radius_max*footprint_radius_max 
	      + instrument.rho_transmitter*instrument.rho_transmitter
	      *(range[i]-instrument.altitude)*(range[i]-instrument.altitude);
	    ms_real eff_footprint_radius2 = footprint_radius2;
	    if (footprint_radius_max2 < footprint_radius2) {
	      eff_footprint_radius2 = footprint_radius_max2;
	    }
	    if (g[i] >= 1.0/(3.0*MS_MU1)) {
	      /* Asymmetry factor is too large for outgoing stream to
		 be scattered directly back to receiver */
	      bscat_inc = transmittance[i]
		*ONE_OVER_FOUR_PI * Iin[it+1][i1] * ssa[i]*ext[i]*2.0
		*(1.0-exp(-(eff_footprint_radius2
			    *Iin[it+1][i1]/Vin[it+1][i1])));
	    }
	    else {
	      bscat_inc = transmittance[i]*ONE_OVER_FOUR_PI
		*(Iin[it+1][i1] * ssa[i]*ext[i]*(1.0+3.0*g[i]*MS_MU1)
		  *(1.0-exp(-(eff_footprint_radius2
			      *Iin[it+1][i1]/Vin[it+1][i1])))
		  +Iout[it+1][i1] * ssa[i]*ext[i]*(1.0-3.0*g[i]*MS_MU1)
		  *(1.0-exp(-(eff_footprint_radius2
			      *Iout[it+1][i1]/Vout[it+1][i1]))));
	    }
	  }
	  else {
	    /* Radar-type Gaussian receiver */
	    if (g[i] >= 1.0/(3.0*MS_MU1)) {
	      bscat_inc = transmittance[i]
		*ONE_OVER_FOUR_PI
		*(Iin[it+1][i1] * ssa[i]*ext[i]*2.0
		  /(1.0+Vin[it+1][i1]/(Iin[it+1][i1]*footprint_radius2)));
	    }
	    else {
	      bscat_inc = transmittance[i]
		*ONE_OVER_FOUR_PI * ssa[i]*ext[i]
		*(Iin[it+1][i1]*(1.0+3.0*g[i]*MS_MU1)
		  /(1.0+Vin[it+1][i1]/(Iin[it+1][i1]*footprint_radius2))
		  +Iout[it+1][i1]*(1.0-3.0*g[i]*MS_MU1)
		  /(1.0+Vout[it+1][i1]/(Iout[it+1][i1]*footprint_radius2)));
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
	int i1 = i+1;
	if (it < n && i == it) {
	  fprintf(stderr, "%d %d %g %g %g %g %g %g %g %g\n", it, i,
		  Iout[it+1][i1], Iin[it+1][i1],
		  Vout[it+1][i1], Vin[it+1][i1],
		  src_power_out[it], src_power_in[it],
		  src_power_out[it]*src_width2[it],
		  src_power_in[it]*src_width2[it]);
	}
	else {
	  fprintf(stderr, "%d %d %g %g %g %g 0 0 0 0\n", it, i,
		  Iout[it+1][i1], Iin[it+1][i1],
		  Vout[it+1][i1], Vin[it+1][i1]);
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
		it, src_power_out[it], it, src_power_in[it], ssa[it], g[it]);
	return MS_FAILURE;
      }
      else if (src_width2[it] < 0.0) {
	fprintf(stderr, "Error: negative source width for two-stream equations at line %d of %s\n",
		__LINE__, __FILE__);
	return MS_FAILURE;
      }
      Iin[it+1][it+1]  += src_power_in[it];
      Iout[it+1][it+1] += src_power_out[it];
      Vin[it+1][it+1]  += src_power_in[it]  * src_width2[it];
      Vout[it+1][it+1] += src_power_out[it] * src_width2[it];
      config->total_src += src_power_in[it] + src_power_out[it];
    }


    config->total_reflected += delta1[0]*Iin[it][1];
  }

  /* Adjoint calculation */
  for (it = nt-1; it >= 0; it--) {
    int max_i = it;
    if (max_i >= n) {
      max_i = n-1;
    }

    /*    max_i = n-1; */
  
    if (it < n) {
      Iin[it+1][it+1]  -= src_power_in[it];
      Iout[it+1][it+1] -= src_power_out[it];
      Vin[it+1][it+1]  -= src_power_in[it]  * src_width2[it];
      Vout[it+1][it+1] -= src_power_out[it] * src_width2[it];

      src_power_in_AD[it] += Iin_next_AD[it+1] 
	+ src_width2[it]*Vin_next_AD[it+1];
      src_power_out_AD[it] += Iout_next_AD[it+1] 
	+ src_width2[it]*Vout_next_AD[it+1];
      src_width2_AD[it] += Vin_next_AD[it+1]*src_power_in[it]
	+ Vout_next_AD[it+1]*src_power_out[it];
    }

    for (i = max_i; i >= 0; i--) {
      ms_real expansion, expansion_AD = 0.0;
      int i1 = i+1;

      for (ifov = 0; ifov < instrument.nfov; ifov++) {
	ms_real footprint_radius2 
	  = instrument.rho_receiver[ifov]*instrument.rho_receiver[ifov]
	  *(range[i]-instrument.altitude)*(range[i]-instrument.altitude);
	int ibscat = ifov*m + (it+i)/2; /* BEST FOR I3RC */
	/* Check that we are not going to overwrite the next point in
	   memory */
	if (ibscat >= (ifov+1)*m) continue;
	if (Vin[it+1][i1] > 0.0 && Vout[it+1][i1] > 0.0) {
	  ms_real bscat_inc_AD = bscat_AD[ibscat];
	  if (config->options & MS_ANNULAR_DETECTORS) {
	    if (ifov < instrument.nfov-1) {
	      bscat_inc_AD -= bscat_AD[ibscat+m];
	    }
	  }
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
	      ms_real exp_term = ONE_OVER_FOUR_PI*2.0
		*exp(-(eff_footprint_radius2*Iin[it+1][i1]/Vin[it+1][i1]));
	      ms_real main_term = ONE_OVER_FOUR_PI*2.0 - exp_term;
	      transmittance_AD[i] += main_term*ext[i]*ssa[i]*Iin[it+1][i1]
		* bscat_inc_AD;
	      ext_AD[i] += main_term*transmittance[i]*ssa[i]*Iin[it+1][i1]
		* bscat_inc_AD;
	      ssa_AD[i] += main_term*transmittance[i]*ext[i]*Iin[it+1][i1]
		* bscat_inc_AD;
	      Iin_next_AD[i1] += ext[i]*ssa[i]*transmittance[i]
		* (main_term + exp_term*eff_footprint_radius2
		   *Iin[it+1][i1]/Vin[it+1][i1])
		* bscat_inc_AD;
	      Vin_next_AD[i1] -= ext[i]*ssa[i]*transmittance[i]
		* exp_term*eff_footprint_radius2*(Iin[it+1][i1]*Iin[it+1][i1])
		* bscat_inc_AD / (Vin[it+1][i1]*Vin[it+1][i1]);

	      /* Note we do not set bscat_AD to zero because original
		 code uses "+=" not "=" */
	    }
	    else {
	      ms_real exp_in = exp(-(eff_footprint_radius2
				     *Iin[it+1][i1]/Vin[it+1][i1]));
	      ms_real exp_out = exp(-(eff_footprint_radius2
				      *Iout[it+1][i1]/Vout[it+1][i1]));
	      
	      transmittance_AD[i] += ONE_OVER_FOUR_PI
		*(Iin[it+1][i1] * ssa[i]*ext[i]*(1.0+3.0*g[i]*MS_MU1)
		  *(1.0-exp_in)
		  +Iout[it+1][i1] * ssa[i]*ext[i]*(1.0-3.0*g[i]*MS_MU1)
		  *(1.0-exp_out)) * bscat_inc_AD;
	      ext_AD[i] += transmittance[i]*ONE_OVER_FOUR_PI
		*(Iin[it+1][i1] * ssa[i]*(1.0+3.0*g[i]*MS_MU1)
		  *(1.0-exp_in)
		  +Iout[it+1][i1] * ssa[i]*(1.0-3.0*g[i]*MS_MU1)
		  *(1.0-exp_out)) * bscat_inc_AD;
	      ssa_AD[i] += transmittance[i]*ONE_OVER_FOUR_PI
		*(Iin[it+1][i1] * ext[i]*(1.0+3.0*g[i]*MS_MU1)
		  *(1.0-exp_in)
		  +Iout[it+1][i1] * ext[i]*(1.0-3.0*g[i]*MS_MU1)
		  *(1.0-exp_out)) * bscat_inc_AD;
	      g_AD[i] += transmittance[i]*ONE_OVER_FOUR_PI
		*(Iin[it+1][i1] * ext[i]*ssa[i]*3.0*MS_MU1
		  *(1.0-exp_in)
		  -Iout[it+1][i1] * ext[i]*ssa[i]*3.0*MS_MU1
		  *(1.0-exp_out)) * bscat_inc_AD;
	      Iin_next_AD[i1] += ext[i]*ssa[i]*transmittance[i]
		* (1.0+3.0*g[i]*MS_MU1) * ONE_OVER_FOUR_PI
		* (1.0-exp_in + exp_in*eff_footprint_radius2
		   *Iin[it+1][i1]/Vin[it+1][i1])
		* bscat_inc_AD;
	      Vin_next_AD[i1] -= ext[i]*ssa[i]*transmittance[i]
		* (1.0+3.0*g[i]*MS_MU1) * ONE_OVER_FOUR_PI
		* exp_in*eff_footprint_radius2
		* Iin[it+1][i1]*Iin[it+1][i1]
		* bscat_inc_AD / (Vin[it+1][i1]*Vin[it+1][i1]);
	      Iout_next_AD[i1] += ext[i]*ssa[i]*transmittance[i]
		* (1.0-3.0*g[i]*MS_MU1) * ONE_OVER_FOUR_PI
		* (1.0-exp_out + exp_out*eff_footprint_radius2
		   *Iout[it+1][i1]/Vout[it+1][i1])
		* bscat_inc_AD;
	      Vout_next_AD[i1] -= ext[i]*ssa[i]*transmittance[i]
		* (1.0-3.0*g[i]*MS_MU1) * ONE_OVER_FOUR_PI
		* exp_out*eff_footprint_radius2
		* Iout[it+1][i1]*Iout[it+1][i1]
		* bscat_inc_AD / (Vout[it+1][i1]*Vout[it+1][i1]);
	    }
	  }
	  else {
	    /* Radar-type Gaussian receiver */
	    if (g[i] >= 1.0/(3.0*MS_MU1)) {
	      ms_real factor = 1.0
		/(1.0+Vin[it+1][i1]/(Iin[it+1][i1]*footprint_radius2));
	      ms_real main_term = ONE_OVER_FOUR_PI*2.0*Iin[it+1][i1]
		* factor;
	      transmittance_AD[i] += main_term * ssa[i]*ext[i]
		*bscat_inc_AD;
	      ext_AD[i] += main_term * ssa[i]*transmittance[i]
		*bscat_inc_AD;
	      ssa_AD[i] += main_term * ext[i]*transmittance[i]
		*bscat_inc_AD;
	      Iin_next_AD[i1] += ONE_OVER_FOUR_PI*2.0
		*transmittance[i]*ext[i]*ssa[i]*bscat_inc_AD
		*(factor + Vin[it+1][i1]*factor*factor
		  /(Iin[it+1][i1]*footprint_radius2));
	      Vin_next_AD[i1] -= ONE_OVER_FOUR_PI*2.0
		*transmittance[i]*ext[i]*ssa[i]*bscat_inc_AD
		*factor*factor/footprint_radius2;
	    }
	    else {
	      ms_real factor_in = 1.0
		/(1.0+Vin[it+1][i1]/(Iin[it+1][i1]*footprint_radius2));
	      ms_real factor_out = 1.0
		/(1.0+Vout[it+1][i1]/(Iout[it+1][i1]*footprint_radius2));
	      ms_real main_in = ONE_OVER_FOUR_PI*Iin[it+1][i1]
		* factor_in;
	      ms_real main_out = ONE_OVER_FOUR_PI*Iout[it+1][i1]
		* factor_out;

	      transmittance_AD[i] += ssa[i]*ext[i]*bscat_inc_AD
		*(main_in*(1.0+3.0*g[i]*MS_MU1)
		  +main_out*(1.0-3.0*g[i]*MS_MU1));
	      ext_AD[i] += ssa[i]*transmittance[i]*bscat_inc_AD
		*(main_in*(1.0+3.0*g[i]*MS_MU1)
		  +main_out*(1.0-3.0*g[i]*MS_MU1));
	      ssa_AD[i] += transmittance[i]*ext[i]*bscat_inc_AD
		*(main_in*(1.0+3.0*g[i]*MS_MU1)
		  +main_out*(1.0-3.0*g[i]*MS_MU1));
	      g_AD[i] += 3.0*MS_MU1*bscat_inc_AD
		*transmittance[i]*ext[i]*ssa[i]
		*(main_in - main_out);
	      Iin_next_AD[i1] += ONE_OVER_FOUR_PI
		*(1.0+3.0*g[i]*MS_MU1)
		*transmittance[i]*ext[i]*ssa[i]*bscat_inc_AD
		*(factor_in + Vin[it+1][i1]*factor_in*factor_in
		  /(Iin[it+1][i1]*footprint_radius2));
	      Vin_next_AD[i1] -= ONE_OVER_FOUR_PI
		*(1.0+3.0*g[i]*MS_MU1)
		*transmittance[i]*ext[i]*ssa[i]*bscat_inc_AD
		*factor_in*factor_in/footprint_radius2;
	      Iout_next_AD[i1] += ONE_OVER_FOUR_PI
		*(1.0-3.0*g[i]*MS_MU1)
		*transmittance[i]*ext[i]*ssa[i]*bscat_inc_AD
		*(factor_out + Vout[it+1][i1]*factor_out*factor_out
		  /(Iout[it+1][i1]*footprint_radius2));
	      Vout_next_AD[i1] -= ONE_OVER_FOUR_PI
		*(1.0-3.0*g[i]*MS_MU1)
		*transmittance[i]*ext[i]*ssa[i]*bscat_inc_AD
		*factor_out*factor_out/footprint_radius2;
	    }
	  }
	}
      }

      expansion = get_expansion_fast(drange, transport_mfp[i],
					     src_width2[i],
					     Iout[it][i1], Vout[it][i1]);
      Vout_AD[i1]       += delta0[i]      * Vout_next_AD[i1];
      delta0_AD[i]      += Vout[it][i1]   * Vout_next_AD[i1];
      Vin_AD[i1]        += delta1[i]      * Vout_next_AD[i1];
      delta1_AD[i]      += Vin[it][i1]    * Vout_next_AD[i1];
      Vout_AD[i1-1]     += delta2[i-1]    * Vout_next_AD[i1];
      delta2_AD[i-1]    += Vout[it][i1-1] * Vout_next_AD[i1];
      Vin_AD[i1+1]      += delta3[i+1]    * Vout_next_AD[i1];
      delta3_AD[i+1]    += Vin[it][i1+1]  * Vout_next_AD[i1];
      Vout_AD[i1+1]     += delta4[i+1]    * Vout_next_AD[i1];
      delta4_AD[i+1]    += Vout[it][i1+1] * Vout_next_AD[i1];
      Vin_AD[i1-1]      += delta5[i-1]    * Vout_next_AD[i1];
      delta5_AD[i-1]    += Vin[it][i1-1]  * Vout_next_AD[i1];
      Iout_AD[i1]       += expansion      * Vout_next_AD[i1];
      expansion_AD =  Iout[it][i1]* Vout_next_AD[i1];
      Vout_next_AD[i1] = 0.0;
      
      get_expansion_fast_AD(drange, transport_mfp[i], src_width2[i],
			    Iout[it][i1], Vout[it][i1],
			    &expansion_AD,
			    transport_mfp_AD+i, src_width2_AD+i,
			    &Iout_AD[i1], &Vout_AD[i1]);
      
      expansion = get_expansion_fast(drange, transport_mfp[i],
				     src_width2[i],
				     Iin[it][i1], Vin[it][i1]);

      Vin_AD[i1]        += delta0[i]      * Vin_next_AD[i1];
      delta0_AD[i]      += Vin[it][i1]    * Vin_next_AD[i1];
      Vout_AD[i1]       += delta1[i]      * Vin_next_AD[i1];
      delta1_AD[i]      += Vout[it][i1]   * Vin_next_AD[i1];
      Vin_AD[i1+1]      += delta2[i+1]    * Vin_next_AD[i1];
      delta2_AD[i+1]    += Vin[it][i1+1]  * Vin_next_AD[i1];
      Vout_AD[i1-1]     += delta3[i-1]    * Vin_next_AD[i1];
      delta3_AD[i-1]    += Vout[it][i1-1] * Vin_next_AD[i1];
      Vin_AD[i1-1]      += delta4[i-1]    * Vin_next_AD[i1];
      delta4_AD[i-1]    += Vin[it][i1-1]  * Vin_next_AD[i1];
      Vout_AD[i1+1]     += delta5[i+1]    * Vin_next_AD[i1];
      delta5_AD[i+1]    += Vout[it][i1+1] * Vin_next_AD[i1];
      Iin_AD[i1]        += expansion      * Vin_next_AD[i1];
      expansion_AD      =  Iin[it][i1]    * Vin_next_AD[i1];
      Vin_next_AD[i1]  = 0.0;
      
      get_expansion_fast_AD(drange, transport_mfp[i], src_width2[i],
			    Iin[it][i1], Vin[it][i1],
			    &expansion_AD,
			    transport_mfp_AD+i, src_width2_AD+i,
			    &Iin_AD[i1], &Vin_AD[i1]);
      
      Iout_AD[i1]       += delta0[i]      * Iout_next_AD[i1];
      delta0_AD[i]      += Iout[it][i1]   * Iout_next_AD[i1];
      Iin_AD[i1]        += delta1[i]      * Iout_next_AD[i1];
      delta1_AD[i]      += Iin[it][i1]    * Iout_next_AD[i1];
      Iout_AD[i1-1]     += delta2[i-1]    * Iout_next_AD[i1];
      delta2_AD[i-1]    += Iout[it][i1-1] * Iout_next_AD[i1];
      Iin_AD[i1+1]      += delta3[i+1]    * Iout_next_AD[i1];
      delta3_AD[i+1]    += Iin[it][i1+1]  * Iout_next_AD[i1];
      Iout_AD[i1+1]     += delta4[i+1]    * Iout_next_AD[i1];
      delta4_AD[i+1]    += Iout[it][i1+1] * Iout_next_AD[i1];
      Iin_AD[i1-1]      += delta5[i-1]    * Iout_next_AD[i1];
      delta5_AD[i-1]    += Iin[it][i1-1]  * Iout_next_AD[i1];

      Iin_AD[i1]        += delta0[i]      * Iin_next_AD[i1];
      delta0_AD[i]      += Iin[it][i1]    * Iin_next_AD[i1];
      Iout_AD[i1]       += delta1[i]      * Iin_next_AD[i1];
      delta1_AD[i]      += Iout[it][i1]   * Iin_next_AD[i1];
      Iin_AD[i1+1]      += delta2[i+1]    * Iin_next_AD[i1];
      delta2_AD[i+1]    += Iin[it][i1+1]  * Iin_next_AD[i1];
      Iout_AD[i1-1]     += delta3[i-1]    * Iin_next_AD[i1];
      delta3_AD[i-1]    += Iout[it][i1-1] * Iin_next_AD[i1];
      Iin_AD[i1-1]      += delta4[i-1]    * Iin_next_AD[i1];
      delta4_AD[i-1]    += Iin[it][i1-1]  * Iin_next_AD[i1];
      Iout_AD[i1+1]     += delta5[i+1]    * Iin_next_AD[i1];
      delta5_AD[i+1]    += Iout[it][i1+1] * Iin_next_AD[i1];
    } /* End of spatial loop */

    for (i = -1; i <= n; i++) {
      int i1 = i+1;
      Iin_next_AD[i1] = Iin_AD[i1];
      Iout_next_AD[i1] = Iout_AD[i1];
      Vin_next_AD[i1] = Vin_AD[i1];
      Vout_next_AD[i1] = Vout_AD[i1];
      Iin_AD[i1] = 0.0;
      Iout_AD[i1] = 0.0;
      Vin_AD[i1] = 0.0;
      Vout_AD[i1] = 0.0;
    }
  } /* End of temporal loop */

  calculate_deltas_AD(n, ext, ssa, g, drange,
		      delta0, delta1, delta2,
		      delta3, delta4, delta5,
		      transmittance, transport_mfp,
		      delta0_AD, delta1_AD, delta2_AD,
		      delta3_AD, delta4_AD, delta5_AD,
			transmittance_AD, transport_mfp_AD,
		      ext_AD, ssa_AD, g_AD);
  
  return MS_SUCCESS;
}
