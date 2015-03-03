/* multiscatter.c -- Adjoint of interface to multiple scattering algorithms

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


   The algorithm is implemented in ANSI C, but a Fortran interface is
   provided.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* For details of how to call the multiple scattering algorithms, see
   the comments in this header file: */
#include "ms.h"
#include "Timer.h"

#include "adept.h"

#ifndef SKIP_NON_ADEPT
#include "cppad/cppad.hpp"
#include "adolc/adolc.h"
#include "Sacado.hpp"
template<> int Sacado::Rad::ADmemblock<double>::n_blocks = 0;


/* This file can access three different automatic differentiation
   libraries by including a header file three times, each time with
   "areal" defined differently */
using CppAD::AD;
#define areal AD<ms_real>
#include "ms_ad.h"
#undef _MS_AD_H
#undef areal

#define areal adouble
#undef _MS_AD_H
#include "ms_ad.h"
#undef areal

#define areal Sacado::Rad::ADvar<double>
#undef _MS_AD_H
#include "ms_ad.h"
#undef areal

#define areal Sacado::ELRFad::DFad<double>
#undef _MS_AD_H
#include "ms_ad.h"
#undef areal
#endif

#define areal adept::aReal
#undef _MS_AD_H
#include "ms_ad.h"

//#define multiscatter simple_algorithm

#ifndef SKIP_NON_ADEPT
typedef Sacado::Rad::ADvar<double> sacado_double;
typedef Sacado::ELRFad::DFad<double> sacado_fad_double;
#endif

/* Perform the multiple scattering calculation, using which ever
   combination of algorithms is appropriate depending on multiple
   scattering regime */
static
int
multiscatter_AD_handcoded(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,    /* Cloud/aerosol equivalent radius, microns */
    const ms_real *ext,       /* Total extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./backscatter ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft., m-1 (NULL for vacuum) */
    const ms_real *ssa_air,   /* Air single-scatter albedo */
    const ms_real *droplet_fraction,     /* Fraction of extinction from droplets */
    const ms_real *pristine_ice_fraction,/* Fraction of ext from pristine ice */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    ms_real *bscat_air_out,   /* Measured backscatter of air, m-1 sr-1 */
    /* Adjoint inputs */
    const ms_real *bscat_AD,
    const ms_real *bscat_air_AD,
    /* Adjoint outputs */
    ms_real *radius_AD,
    ms_real *ext_AD,
    ms_real *ssa_AD,
    ms_real *g_AD,
    ms_real *ext_bscat_ratio_AD,
    ms_real *ext_air_AD,
    ms_real *droplet_fraction_AD,
    ms_real *pristine_ice_fraction_AD)
{
  ms_instrument instrument_1fov = instrument;
  ms_real ss_multiplier = 1.0;
  int status = MS_SUCCESS;
  int i;
  int ifov;
  int require_bscat_air = (bscat_air_out != NULL);
  int nan_input = 0;

  if (config->options & MS_CHECK_FOR_NAN
      && ms_isnan(n, ext_AD)) {
    fprintf(stderr, "Warning: ext_AD input to multiscatter_AD() contains NaNs\n");
    nan_input = 1;
  }
  //  fprintf(stderr, "!!!");
    for (i = 0; i < n; i++) {
      //      fprintf(stderr, " %g", ext_AD[i]);
      ext_AD[i] = 0.0;
    }
    //  fprintf(stderr, "\n");

  /* This variable contains the same instrument settings, except for
     having only one field-of-view; it is modified when we loop over
     the fields-of-view */
  instrument_1fov.nfov = 1;

  /* Work out single-scattering factor: this ensures that the first
     field-of-view is properly calibrated and that the backscatter in
     the subsequent fields-of-view are calibrated to the first such
     that the backscatter is proportional to the energy received. */
  if (instrument.receiver_type == MS_TOP_HAT) {
    ss_multiplier
      = 1.0/(1.0-exp(-instrument.rho_receiver[0]*instrument.rho_receiver[0]/
		     (instrument.rho_transmitter*instrument.rho_transmitter)));
  }
  else {
    ss_multiplier
      = 1.0 + instrument.rho_transmitter*instrument.rho_transmitter
      /(instrument.rho_receiver[0]*instrument.rho_receiver[0]);
  }
  config->ss_multiplier = ss_multiplier;

  if (!(config->options & MS_QUIET) && ss_multiplier != 1.0) {
    fprintf(stderr, "Factor to ensure calibration of first FOV: %g\n",
	    ss_multiplier);
  }


  /* FIRST DO SMALL-ANGLE SCATTERING CALCULATION */

  if (config->small_angle_algorithm != MS_SINGLE_AND_SMALL_ANGLE_NONE) {
    /* In case of annular detectors we may need to manipulate the
       adjoints before passing them into the functions - here we
       allocate the storage for the manipulated vectors */
    ms_real bscat_annular_AD_data[m];
    ms_real bscat_air_annular_AD_data[m];

    /* For small-angle scattering we loop over the fields-of-view */
    for (ifov = 0; ifov < instrument.nfov; ifov++) {
      /* Set the receiver field of view appropriately */
      const ms_real* bscat_annular_AD = bscat_AD + ifov*m;
      const ms_real* bscat_air_annular_AD 
	= bscat_air_AD + require_bscat_air*ifov*m;
      ms_real rho_receiver = instrument.rho_receiver[ifov];
      ms_real fov_ss_multiplier = 1;

      instrument_1fov.rho_receiver = &rho_receiver;

      if (ifov > 0) {
	if (instrument.receiver_type == MS_TOP_HAT) {
	  fov_ss_multiplier
	    = ss_multiplier*(1.0-exp(-instrument.rho_receiver[ifov]
				     *instrument.rho_receiver[ifov]/
				     (instrument.rho_transmitter
				      *instrument.rho_transmitter)));
	}
	else {
	  fov_ss_multiplier = ss_multiplier
	    / (1.0 + instrument.rho_transmitter*instrument.rho_transmitter
	       /(instrument.rho_receiver[ifov]*instrument.rho_receiver[ifov]));
	}
      }


      /* We need to manipulate adjoints if we are using annular
	 detectors AND we have more than one field of view AND it is
	 not the final field of view */
      if (config->options & MS_ANNULAR_DETECTORS) {
	if (ifov < instrument.nfov-1) {
	  int i;
	  for (i = 0; i < m; i++) {
	    bscat_annular_AD_data[i] = 1.0/*fov_ss_multiplier*/
	      * (bscat_annular_AD[i] - bscat_annular_AD[i+m]);
	  }
	  bscat_annular_AD = bscat_annular_AD_data;
	  if (require_bscat_air) {
	    for (i = 0; i < m; i++) {
	      bscat_air_annular_AD_data[i] = 1.0/*fov_ss_multiplier*/
		* (bscat_air_annular_AD[i] - bscat_air_annular_AD[i+m]);
	    }
	    bscat_air_annular_AD = bscat_air_annular_AD_data;
	  }
	}
	/*
	else {
	  int i;
	  for (i = 0; i < m; i++) {
	    bscat_annular_AD_data[i] = fov_ss_multiplier
	      * bscat_annular_AD[i];
	  }
	  bscat_annular_AD = bscat_annular_AD_data;
	  if (require_bscat_air) {
	    for (i = 0; i < m; i++) {
	      bscat_air_annular_AD_data[i] = fov_ss_multiplier
		* bscat_air_annular_AD[i];
	    }
	    bscat_air_annular_AD = bscat_air_annular_AD_data;
	  }
	}
	*/

      }

      else if (fov_ss_multiplier != 1.0
	       && ifov > 0) {
	int i;
	for (i = 0; i < m; i++) {
	  bscat_annular_AD_data[i] = fov_ss_multiplier
	    * bscat_annular_AD[i];
	}
	bscat_annular_AD = bscat_annular_AD_data;
	if (require_bscat_air) {
	  for (i = 0; i < m; i++) {
	    bscat_air_annular_AD_data[i] = fov_ss_multiplier
	      * bscat_air_annular_AD[i];
	  }
	  bscat_air_annular_AD = bscat_air_annular_AD_data;
	}
      }

      if (config->small_angle_algorithm == MS_SINGLE_SCATTERING) {
	/* If no forward lobe (radar case) then we have no small-angle
	   multiple scattering contribution and we simply calculate
	   the single scattering here, with wide-angle scattering
	   later. */
	status
	  = ms_singlescatter_AD(n, instrument_1fov, surface,
				range, ext, ext_bscat_ratio,
				ext_air,
				bscat_out + ifov*m,
				bscat_air_out + ifov*m*require_bscat_air,
				bscat_annular_AD,
				bscat_air_annular_AD,
				ext_AD, ext_bscat_ratio_AD);
      }
      else {
	if (config->small_angle_algorithm != MS_SMALL_ANGLE_PVC_EXPLICIT) {
	  /* Perform small-angle multiple-scattering calculation using an
	     efficient algorithm */
	  if (config->small_angle_algorithm == MS_SMALL_ANGLE_PVC_FAST_LAG) {
	    fprintf(stderr, "Error: adjoint not available with small-angle lag calculation\n");
	    return MS_FAILURE;
	  }
	  else if (config->small_angle_algorithm == MS_SMALL_ANGLE_PVC_FAST) {
	    status = ms_fast_small_angle_AD(n, 
			 instrument_1fov, surface, 
			 range, radius, ext, ext_bscat_ratio, ext_air,
			 droplet_fraction, pristine_ice_fraction,
			 bscat_out + ifov*m,
			 bscat_air_out + ifov*m*require_bscat_air,
		         bscat_AD, bscat_air_AD, radius_AD, 
					    ext_AD, ext_bscat_ratio_AD, NULL, /* Don't do ext_air_AD */
			 droplet_fraction_AD, pristine_ice_fraction_AD);
	  }
	  else {
	    fprintf(stderr, "Error: adjoint not available for O(N^2) small-angle calculation (Hogan 2006): use the \"-fast-sa\" option\n");
	    return MS_FAILURE;
	  }
	}
	else {
	  fprintf(stderr, "Error: adjoint not available for explicit small-angle calculation: use the \"-fast-sa\" option\n");
	  return MS_FAILURE;
	}
	if (status != MS_SUCCESS) {
	  return status;
	}

      } /* End of if single scattering else small-angle scattering */

      /* To ensure the correct calibration, we may need to scale the
	 wider fields of view - here we multiply by the multiplier for
	 the first field-of-view, and divide by the multiplier
	 calculated for the current field-of-view.  At the wide-angle
	 scattering stage, we simply multiply by the multiplier for
	 the first field-of-view.  This ensures that for a given
	 field-of-view, the narrow-angle and wide-angle contributions
	 are correctly scaled relative to each other (normally one
	 would only multiply the wide-angle component by the
	 multiplier for the current FOV), and relative to the first
	 field-of-view */
      if (ifov > 0) {
	ms_real* bscat_cur = bscat_out+ifov*m;
	ms_real* bscat_air_cur = bscat_air_out+require_bscat_air*ifov*m;
	for (i = 0; i < n; i++, bscat_cur++) {
	  (*bscat_cur) *= fov_ss_multiplier;
	}
	if (bscat_air_cur) {
	  for (i = 0; i < n; i++, bscat_air_cur++) {
	    (*bscat_air_cur) *= fov_ss_multiplier;
	  }
	}
	
	if (!(config->options & MS_QUIET) && ss_multiplier != 1.0) {
	  fprintf(stderr, "   Equivalent factor for FOV %d: %g\n",
		  ifov+1, ss_multiplier/fov_ss_multiplier);
	}
      }
    } /* Loop over fields-of-view */

    if (config->options & MS_ANNULAR_DETECTORS) {
      /* Need to subtract backscatter contributions from each other */
      for (ifov = instrument.nfov-1; ifov > 0; ifov--) {
	int i;
	int i_inner = (ifov-1)*m;
	int i_outer = ifov*m;
	for (i = 0; i < m; i++, i_inner++, i_outer++) {
	  bscat_out[i_outer] -= bscat_out[i_inner];
	}
	if (bscat_air_out) {
	  i_inner = (ifov-1)*m;
	  i_outer = ifov*m;
	  for (i = 0; i < m; i++, i_inner++, i_outer++) {
	    bscat_air_out[i_outer] -= bscat_air_out[i_inner];
	  }
	}
      }
    }
  } /* If perform small-angle calculation */

  /* Return if status is non-zero (indicating that an error
     occurred) or if ssa or g are not set (indicating that the
     wide-angle calculation is not to be performed). */
  if (status != MS_SUCCESS
      || config->wide_angle_algorithm == MS_WIDE_ANGLE_NONE) {
    return status;
  }
  if ( !ssa || !g ) {
    fprintf(stderr, "Error: tried to call wide-angle part of calculation without single scattering albedo or asymmetry factor being set\n");
    return MS_FAILURE;
  }

  if (config->small_angle_algorithm == MS_SINGLE_AND_SMALL_ANGLE_NONE) {
    /* Reset backscatter to zero, since wide-angle calculation
       increments existing values */
    int i;
    for (i = 0; i < m*instrument.nfov; i++) {
      bscat_out[i] = 0.0;
    }
  }

  /* NOW DO WIDE-ANGLE SCATTERING CALCULATION WITH ADJOINT */
  status = ms_wide_angle_AD(n, m,
	    config, instrument, surface, range, radius,
	    ext, ssa, g, ext_air, ssa_air,
	    bscat_out,
	    bscat_AD, 
            radius_AD, ext_AD, ssa_AD, g_AD);

  if (status != MS_SUCCESS) {
    return status;
  }

  if (config->options & MS_CHECK_FOR_NAN
      && !nan_input
      && ms_isnan(n, ext_AD)) {
    fprintf(stderr, "Warning: algorithm produced NaNs in ext_AD\n");
    ms_dump_input(NULL, n, config, instrument, surface, range, radius,
		  ext, ssa, g, ext_bscat_ratio, ext_air, ssa_air,
		  droplet_fraction, pristine_ice_fraction,
		  bscat_AD, bscat_air_AD);
  }
  return MS_SUCCESS;
}
Timer timer(3);

int
multiscatter_AD(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,    /* Cloud/aerosol equivalent radius, microns */
    const ms_real *ext,       /* Total extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./backscatter ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft., m-1 (NULL for vacuum) */
    const ms_real *ssa_air,   /* Air single-scatter albedo */
    const ms_real *droplet_fraction,     /* Fraction of extinction from droplets */
    const ms_real *pristine_ice_fraction,/* Fraction of ext from pristine ice */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    ms_real *bscat_air_out,   /* Measured backscatter of air, m-1 sr-1 */
    /* Adjoint inputs */
    const ms_real *bscat_AD,
    const ms_real *bscat_air_AD,
    /* Adjoint outputs */
    ms_real *radius_AD,
    ms_real *ext_AD,
    ms_real *ssa_AD,
    ms_real *g_AD,
    ms_real *ext_bscat_ratio_AD,
    ms_real *ext_air_AD,
    ms_real *droplet_fraction_AD,
    ms_real *pristine_ice_fraction_AD) {
  static bool is_printed = false;
  //  static bool jac_is_printed = false;
  static bool jac_is_printed = true;

  if (config->options & MS_AUTOMATIC_DIFFERENTIATION_ADEPT) {
    static adept::Stack ad_stack;
    ad_stack.preallocate_statements(1000000);
    ad_stack.preallocate_operations(4000000);


    areal a_radius[n];
    areal a_ext[n];
    areal a_ssa[n];
    areal a_g[n];
    areal a_ext_bscat_ratio[n];
    areal a_droplet_fraction[n];
    areal a_pristine_ice_fraction[n];

    areal a_bscat_out[n];
    areal a_bscat_air_tmp[n];
    areal* a_bscat_air_out = 0;
    if (bscat_air_out) {
      a_bscat_air_out = a_bscat_air_tmp;
    }

    set_values(a_radius, n, radius);
    set_values(a_ext, n, ext);
    set_values(a_ssa, n, ssa);
    set_values(a_g, n, g);
    set_values(a_ext_bscat_ratio, n, ext_bscat_ratio);
    set_values(a_droplet_fraction, n, droplet_fraction);
    set_values(a_pristine_ice_fraction, n, pristine_ice_fraction);

    ad_stack.new_recording();

    timer.start(0);
    int status = multiscatter(n, m,config, instrument, surface, range,
			      a_radius, a_ext, a_ssa, a_g, a_ext_bscat_ratio,
			      ext_air, ssa_air,
			      a_droplet_fraction, a_pristine_ice_fraction,
			      a_bscat_out, a_bscat_air_out);
    timer.start(1);
    set_gradients(a_bscat_out, n, bscat_AD);
    if (bscat_air_out && bscat_air_AD) {
      set_gradients(a_bscat_air_out, n, bscat_air_AD);
    }

    if (config->options & MS_JACOBIAN) {
      ad_stack.independent(a_radius, n);
      ad_stack.independent(a_ext, n);
      ad_stack.independent(a_ssa, n);
      ad_stack.independent(a_g, n);
      ad_stack.independent(a_ext_bscat_ratio, n);
      ad_stack.independent(a_droplet_fraction, n);
      ad_stack.independent(a_pristine_ice_fraction, n);
      ad_stack.dependent(a_bscat_out, n);
      if (bscat_air_out && bscat_air_AD) {
	ad_stack.dependent(a_bscat_air_out, n);
      }
      int nout = 50*(1+(bscat_air_out != 0));
      ms_real jacobian[nout*7*n];
      if (config->force_jacobian > 0) {
	ad_stack.jacobian_forward(jacobian);
      }
      else if (config->force_jacobian < 0) {
	ad_stack.jacobian_reverse(jacobian);
      }
      else {
	ad_stack.jacobian(jacobian);
      }
      if (!jac_is_printed) {
	for (int i = 0; i < ad_stack.n_independent(); i++) {
	  for (int j = 0; j < ad_stack.n_dependent(); j++) {
	    std::cerr << " " << jacobian[i*ad_stack.n_dependent() + j];
	  }
	  std::cerr << "\n";
	}
	jac_is_printed = true;
      }
    }
    else {
      ad_stack.compute_adjoint();

      get_gradients(a_radius, n, radius_AD);
      get_gradients(a_ext, n, ext_AD);
      get_gradients(a_ssa, n, ssa_AD);
      get_gradients(a_g, n, g_AD);
      get_gradients(a_ext_bscat_ratio, n, ext_bscat_ratio_AD);
      get_gradients(a_droplet_fraction, n, droplet_fraction_AD);
      get_gradients(a_pristine_ice_fraction, n, pristine_ice_fraction_AD);
    
      get_values(a_bscat_out, n, bscat_out);
      if (bscat_air_out) {
	get_values(a_bscat_air_out, n, bscat_air_out);
      }
    }

    timer.stop();

    if (!is_printed) {
      std::cerr << ad_stack;
      is_printed = true;
    }
    
    return status;

  }
#ifndef SKIP_NON_ADEPT
  else if (config->options & MS_AUTOMATIC_DIFFERENTIATION_CPPAD) {
    if (!is_printed) {
      CppAD::thread_alloc::hold_memory(true);
    }

    std::vector<AD<ms_real> > X(7*n);
    for (int i = 0; i < n; i++) {
      X[i] = radius[i];
      X[i+n] = ext[i];
      X[i+2*n] = ssa[i];
      X[i+3*n] = g[i];
      X[i+4*n] = ext_bscat_ratio[i];
      X[i+5*n] = droplet_fraction[i];
      X[i+6*n] = pristine_ice_fraction[i];
    }
    CppAD::Independent(X);

    int nout = 2*n;
    if (!bscat_air_out) {
      nout = n;
    }
    std::vector<AD<ms_real> > Y(nout);
    AD<ms_real>* bscat_air_out_active = 0;
    if (bscat_air_out) {
      bscat_air_out_active = &Y[n];
    }
    timer.start(0);
    int status = multiscatter(n,m,config,instrument,surface, range,
			      &X[0], &X[n], &X[2*n], &X[3*n], &X[4*n],
			      ext_air, ssa_air,
			      &X[5*n], &X[6*n],
			      &Y[0], bscat_air_out_active);
    timer.start(1);
    if (!is_printed) {
      
      // I needed to hack the CppAD library for the following line to work:
      size_t mem = X[0].tape_this()->memory();
      std::cerr << "Memory used by tape: " << mem << " bytes\n";
    }

    CppAD::ADFun<ms_real> f(X, Y);
    if (!is_printed) {
      // I needed to hack the CppAD library for the following line to work:
      //      size_t mem = Y[10].tape_this()->memory();
      //      std::cerr << "Memory used by tape: " << mem << " bytes\n";
      
      std::cerr << "Memory used by function: " << f.size_op_seq() << " bytes\n";
    }

    std::vector<ms_real> w(nout), dw(7*n);
    for (int i = 0; i < n; i++) {
      w[i] = bscat_AD[i];
      if (bscat_air_out && bscat_air_AD) {
	w[i+n] = bscat_air_AD[i];
      }
    }
    timer.start(2);
    if (config->options & MS_JACOBIAN) {
      static std::vector<double> jacobian(nout * 7*n);
      //      static std::vector<double> jacobian;
      static std::vector<double> x(7*n);
      for (int i = 0; i < n; i++) {
	x[i] = radius[i];
	x[i+n] = ext[i];
	x[i+2*n] = ssa[i];
	x[i+3*n] = g[i];
	x[i+4*n] = ext_bscat_ratio[i];
	x[i+5*n] = droplet_fraction[i];
	x[i+6*n] = pristine_ice_fraction[i];
      }

      if (config->force_jacobian > 0) {
	CppAD::JacobianFor(f, x, jacobian);
      }
      else if (config->force_jacobian < 0) {
	CppAD::JacobianRev(f, x, jacobian);
      }
      else {
	jacobian = f.Jacobian(x);
      }
      if (!jac_is_printed) {
	for (int i = 0; i < 7*n; i++) {
	  for (int j = 0; j < nout; j++) {
	    std::cerr << " " << jacobian[i + j*7*n];
	  }
	  std::cerr << "\n";
	}
	jac_is_printed = true;
      }

    }
    else {
      dw = f.Reverse(1, w);
    
      for (int i = 0; i < n; i++) {
	radius_AD[i] = dw[i];
	ext_AD[i] = dw[i+n];
	ssa_AD[i] = dw[i+2*n];
	g_AD[i] = dw[i+3*n];
	ext_bscat_ratio_AD[i] = dw[i+4*n];
	droplet_fraction_AD[i] = dw[i+5*n];
	pristine_ice_fraction_AD[i] = dw[i+6*n];      
      }
    }
    timer.stop();

    is_printed = true;

    return status;
    
  }
  else if (config->options & MS_AUTOMATIC_DIFFERENTIATION_ADOLC) {
    timer.start(0);
    adouble a_radius[n];
    adouble a_ext[n];
    adouble a_ssa[n];
    adouble a_g[n];
    adouble a_ext_bscat_ratio[n];
    adouble a_droplet_fraction[n];
    adouble a_pristine_ice_fraction[n];

    adouble a_bscat_out[n];
    adouble a_bscat_air_tmp[n];
    adouble* a_bscat_air_out = 0;
    if (bscat_air_out) {
      a_bscat_air_out = a_bscat_air_tmp;
    }
    trace_on(1,1);
    for (int i = 0; i < n; i++) {
      a_radius[i] <<= radius[i];
    }
    for (int i = 0; i < n; i++) {
      a_ext[i] <<= ext[i];
    }
    for (int i = 0; i < n; i++) {
      a_ssa[i] <<= ssa[i];
    }
    for (int i = 0; i < n; i++) {
      a_g[i] <<= g[i];
    }
    for (int i = 0; i < n; i++) {
      a_ext_bscat_ratio[i] <<= ext_bscat_ratio[i];
    }
    for (int i = 0; i < n; i++) {
      a_droplet_fraction[i] <<= droplet_fraction[i];
    }
    for (int i = 0; i < n; i++) {
      a_pristine_ice_fraction[i] <<= pristine_ice_fraction[i];
    }
    int status = multiscatter(n, m,config, instrument, surface, range,
			      a_radius, a_ext, a_ssa, a_g, a_ext_bscat_ratio,
			      ext_air, ssa_air,
			      a_droplet_fraction, a_pristine_ice_fraction,
			      a_bscat_out, a_bscat_air_out);
    for (int i = 0; i < n; i++) {
      a_bscat_out[i] >>= bscat_out[i];
      if (bscat_air_out && bscat_air_AD) {
	a_bscat_air_out[i] >>= bscat_air_out[i];
      }
    }
    trace_off();
    int nout = n;
    if (bscat_air_out && bscat_air_AD) {
      nout *= 2;
    }
    ms_real w[nout];

    for (int i = 0; i < n; i++) {
      w[i] = bscat_AD[i];
      if (bscat_air_out && bscat_air_AD) {
	w[i+n] = bscat_air_AD[i];
      }
    }

    if (!is_printed) {
      size_t counts[11];
      tapestats(1, counts);
      std::cerr << "Counts: " 
		<< counts[0] << " "
		<< counts[1] << " "
		<< counts[2] << " "
		<< counts[3] << " "
		<< counts[4] << " "
		<< counts[5] << "\n";
    }
    timer.start(2);

    if (config->options & MS_JACOBIAN) {
      static ms_real** jac = myalloc2(nout, 7*n);
      ms_real x[7*n];
      for (int i = 0; i < n; i++) {
	x[i] = radius[i];
	x[i+n] = ext[i];
	x[i+2*n] = ssa[i];
	x[i+3*n] = g[i];
	x[i+4*n] = ext_bscat_ratio[i];
	x[i+5*n] = droplet_fraction[i];
	x[i+6*n] = pristine_ice_fraction[i];
      }
      if (config->force_jacobian < 0) {
	double** I = myallocI2(nout);
	double* result = myalloc1(nout);
	int rc = zos_forward(1, nout, 7*n, 1, x, result);
	if (rc < 0) {
	  std::cerr << "ERROR OCCURRED IN ADOL-C's zof_forward()\n";
	  exit(rc);
	}
	MINDEC(rc,fov_reverse(1, nout, 7*n, nout, I, jac));
	myfreeI2(nout, I);
	myfree1(result);
      }
      else if (config->force_jacobian > 0) {
	double* result = myalloc1(nout);
	double** I = myallocI2(7*n);
	int rc = fov_forward(1, nout, 7*n, 7*n, x, I, result, jac);
	myfreeI2(7*n, I);
	myfree1(result);
      }
      else {
	jacobian(1, nout, 7*n, x, jac);
      }
      if (!jac_is_printed) {
	for (int i = 0; i < 7*n; i++) {
	  for (int j = 0; j < nout; j++) {
	    std::cerr << " " << jac[j][i];
	  }
	  std::cerr << "\n";
	}
	jac_is_printed = true;
      }

    }
    else {
      ms_real  dw[7*n];
      reverse(1,nout,7*n,0,w,dw);        /* evaluate the (i+1)-st deriv. */
      timer.stop();
      for (int i = 0; i < n; i++) {
	radius_AD[i] = dw[i];
	ext_AD[i] = dw[i+n];
	ssa_AD[i] = dw[i+2*n];
	g_AD[i] = dw[i+3*n];
	ext_bscat_ratio_AD[i] = dw[i+4*n];
	droplet_fraction_AD[i] = dw[i+5*n];
	pristine_ice_fraction_AD[i] = dw[i+6*n];      
      }
    }

    is_printed = true;
    return status;
  }


  else if (config->options & MS_AUTOMATIC_DIFFERENTIATION_SACADO) {
    sacado_double a_radius[n];
    sacado_double a_ext[n];
    sacado_double a_ssa[n];
    sacado_double a_g[n];
    sacado_double a_ext_bscat_ratio[n];
    sacado_double a_droplet_fraction[n];
    sacado_double a_pristine_ice_fraction[n];

    sacado_double a_bscat_out[n];
    sacado_double a_bscat_air_tmp[n];
    sacado_double* a_bscat_air_out = 0;

    if (bscat_air_out) {
      a_bscat_air_out = a_bscat_air_tmp;
    }

    int nout = 2*n;
    if (!bscat_air_out) {
      nout = n;
    }
    /*
    std::vector<sacado_double> Y(nout);
    sacado_double* bscat_air_out_active = 0;
    if (bscat_air_out) {
      bscat_air_out_active = &Y[n];
    }
    sacado_double bscat_out_active = &Y[0];

    }
    */

    for (int i = 0; i < n; i++) {
      a_radius[i] = radius[i];
      a_ext[i] = ext[i];
      a_ssa[i] = ssa[i];
      a_g[i] = g[i];
      a_ext_bscat_ratio[i] = ext_bscat_ratio[i];
      a_droplet_fraction[i] = droplet_fraction[i];
      a_pristine_ice_fraction[i] = pristine_ice_fraction[i];
    }

    timer.start(0);
    int status = multiscatter(n, m,config, instrument, surface, range,
			      a_radius, a_ext, a_ssa, a_g, a_ext_bscat_ratio,
			      ext_air, ssa_air,
			      a_droplet_fraction, a_pristine_ice_fraction,
			      a_bscat_out, a_bscat_air_out);
    timer.start(1);


    if (config->options & MS_JACOBIAN) {
      std::cerr << "NO SACADO JACOBIAN YET\n";
      exit(1);
      /*
      ad_stack.independent(a_radius, n);
      ad_stack.independent(a_ext, n);
      ad_stack.independent(a_ssa, n);
      ad_stack.independent(a_g, n);
      ad_stack.independent(a_ext_bscat_ratio, n);
      ad_stack.independent(a_droplet_fraction, n);
      ad_stack.independent(a_pristine_ice_fraction, n);
      ad_stack.dependent(a_bscat_out, n);
      if (bscat_air_out && bscat_air_AD) {
	ad_stack.dependent(a_bscat_air_out, n);
      }
      int nout = 50*(1+(bscat_air_out != 0));
      ms_real jacobian[nout*7*n];
      ad_stack.jacobian(jacobian);
      if (!jac_is_printed) {
	for (int i = 0; i < ad_stack.n_independent(); i++) {
	  for (int j = 0; j < ad_stack.n_dependent(); j++) {
	    std::cerr << " " << jacobian[i*ad_stack.n_dependent() + j];
	  }
	  std::cerr << "\n";
	}
	jac_is_printed = true;
      }
      */
    }
    else {
    
    sacado_double objective_func = 0.0;
    for (int i = 0; i < n; i++) {
      objective_func += a_bscat_out[i] * bscat_AD[i];
    }
    if (bscat_air_out && bscat_air_AD) {
      for (int i = 0; i < n; i++) {
	objective_func += a_bscat_air_out[i] * bscat_air_AD[i];
      }
    }
      Sacado::Rad::ADvar<double>::Gradcomp();
    
      /*
      sacado_double* a_y[nout];
      double y_AD[nout];
      for (int i = 0; i < n; i++) {
	a_y[i] = &a_bscat_out[i];
	y_AD[i] = bscat_AD[i];
      }
      if (bscat_air_out) {
	for (int i = 0; i < n; i++) {
	  a_y[i+n] = &a_bscat_air_out[i];
	  y_AD[i+n] = bscat_air_AD[i];
	}
      }
      Sacado::Rad::ADvar<double>::Weighted_Gradcomp(nout, a_y, y_AD);
      */
      for (int i = 0; i < n; i++) {
	radius_AD[i] = a_radius[i].adj();
	ext_AD[i] = a_ext[i].adj();
	ssa_AD[i] = a_ssa[i].adj();
	g_AD[i] = a_g[i].adj();
	ext_bscat_ratio_AD[i] = a_ext_bscat_ratio[i].adj();
	droplet_fraction_AD[i] = a_droplet_fraction[i].adj();
	pristine_ice_fraction_AD[i] = a_pristine_ice_fraction[i].adj();
	bscat_out[i] = a_bscat_out[i].val();
	if (bscat_air_out) {
	  bscat_air_out[i] = a_bscat_air_out[i].val();
	}
      }
      
      timer.stop();
      
      /*
	if (!is_printed) {
	std::cerr << ad_stack;
	is_printed = true;
	}
      */
   
     }
    return status;  
 }


  else if (config->options & MS_AUTOMATIC_DIFFERENTIATION_SACADO_FAD) {
    sacado_fad_double a_radius[n];
    sacado_fad_double a_ext[n];
    sacado_fad_double a_ssa[n];
    sacado_fad_double a_g[n];
    sacado_fad_double a_ext_bscat_ratio[n];
    sacado_fad_double a_droplet_fraction[n];
    sacado_fad_double a_pristine_ice_fraction[n];

    sacado_fad_double a_bscat_out[n];
    sacado_fad_double a_bscat_air_tmp[n];
    sacado_fad_double* a_bscat_air_out = 0;

    if (bscat_air_out) {
      a_bscat_air_out = a_bscat_air_tmp;
    }

    int nout = 2*n;
    if (!bscat_air_out) {
      nout = n;
    }
    /*
    std::vector<sacado_fad_double> Y(nout);
    sacado_fad_double* bscat_air_out_active = 0;
    if (bscat_air_out) {
      bscat_air_out_active = &Y[n];
    }
    sacado_fad_double bscat_out_active = &Y[0];

    }
    */

    int nx = 7*n;

    for (int i = 0; i < n; i++) {
      a_radius[i] = radius[i];
      a_ext[i] = ext[i];
      a_ssa[i] = ssa[i];
      a_g[i] = g[i];
      a_ext_bscat_ratio[i] = ext_bscat_ratio[i];
      a_droplet_fraction[i] = droplet_fraction[i];
      a_pristine_ice_fraction[i] = pristine_ice_fraction[i];
    }
    for (int i = 0; i < n; i++) {
      a_radius[i].resize(nx);
      a_radius[i].fastAccessDx(i) = 1.0;
      a_ext[i].resize(nx);
      a_ext[i].fastAccessDx(i+n) = 1.0;
      a_ssa[i].resize(nx);
      a_ssa[i].fastAccessDx(i+2*n) = 1.0;
      a_g[i].resize(nx);
      a_g[i].fastAccessDx(i+3*n) = 1.0;
      a_ext_bscat_ratio[i].resize(nx);
      a_ext_bscat_ratio[i].fastAccessDx(i+4*n) = 1.0;
      a_droplet_fraction[i].resize(nx);
      a_droplet_fraction[i].fastAccessDx(i+5*n) = 1.0;
      a_pristine_ice_fraction[i].resize(nx);
      a_pristine_ice_fraction[i].fastAccessDx(i+6*n) = 1.0;
    }


    timer.start(0);
    int status = multiscatter(n, m,config, instrument, surface, range,
			      a_radius, a_ext, a_ssa, a_g, a_ext_bscat_ratio,
			      ext_air, ssa_air,
			      a_droplet_fraction, a_pristine_ice_fraction,
			      a_bscat_out, a_bscat_air_out);
    timer.stop();

      if (!jac_is_printed) {
	for (int i = 0; i < n; i++) {
	  std::cerr << "? ";
	  for (int j = 0; j < n; j++) {
	    std::cerr << " " << a_bscat_out[j].dx(i+n);
	  }
	  std::cerr << "\n";
	}
	jac_is_printed = true;
      }
    return status;  

  }
#endif


  else {
    timer.start(0);
    int status = multiscatter_AD_handcoded(n, m, config, instrument, surface, range, radius,
					   ext, ssa, g, ext_bscat_ratio, ext_air, ssa_air,
					   droplet_fraction, pristine_ice_fraction,
					   bscat_out, bscat_air_out,
					   bscat_AD, bscat_air_AD,
					   radius_AD, ext_AD, ssa_AD, g_AD, 
					   ext_bscat_ratio_AD, ext_air_AD, 
					   droplet_fraction_AD, pristine_ice_fraction_AD);
      timer.stop();
      return status;
  }
}
