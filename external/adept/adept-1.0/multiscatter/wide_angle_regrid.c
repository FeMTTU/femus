/* wide_angle_regrid.c -- Regrid before wide-angle multiple scattering

   Copyright (C) 2011 Robin Hogan <r.j.hogan@reading.ac.uk> 

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

/* Do the wide-angle part of the calculation after regridding on to a
   regular grid (if necessary) or skipping gates in which wide-angle
   multiple scattering is likely to be negligible */
int ms_wide_angle_regrid(
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
    ms_real *bscat_out)       /* Measured backscatter, m-1 sr-1 */
{
  ms_real optical_depth_to_start_2way = 0.0;
  int first_gate = config->first_wide_angle_gate;
  if (first_gate < 0) {
    /* Compare mean-free-path to footprint to work out when to start
       wide-angle calculation */
    int i;
    for (i = 0; i < n-1; i++) {
      if (ext[i] > 0.0 && ssa[i] > 0.0) {
	/* The footprint radius at the altitude of the cloud of the
	   widest field of view */
	ms_real footprint = fabs(instrument.altitude-range[i])
	  *instrument.rho_receiver[instrument.nfov-1];
	/* The relevant mean-free-path */
	ms_real mfp = 1.0/(ext[i]*ssa[i]*(1.0 - g[i]));
	if (mfp < footprint*MS_MFP_FOOTPRINT_THRESHOLD) {
	  /* Threshold met: start multiple scattering calculation from
	     this gate */
	  first_gate = i;
	  break;
	}
      }
      /* Estimate optical depth to the start of the multiple scattering
	 region: this will be used in calculating the source terms for the
	 TDTS calculation later on */
      ms_real drange = ms_get_drange(i, n, range);
      if (ext_air) {
	optical_depth_to_start_2way += 2.0*ext_air[i]*drange;
      }
      if (config->wide_angle_algorithm == MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE) {
	/* If there is no forward lobe then we simply add the two-way
	   optical depth of the layer to the total */
	optical_depth_to_start_2way += 2.0*ext[i]*drange;
      }
      else {
	/* There is a forward lobe - this means we need to replicate
	   the calculations in wide_angle.c: use "diffraction scaling"
	   on the outward journey and "delta-Eddington scaling" on the
	   inward journey; odd but it works... */
	/* Outward journey: assume half the extinguished energy is
	   forward scattered and behaves as if it had not been
	   scattered at all */
	ms_real od_layer_outward = 0.5*ext[i]*drange;
	/* Inward journey: use traditional delta-Eddington */
	ms_real od_layer_inward = ext[i]*(1.0-ssa[i]*g[i]*g[i])*drange;
	optical_depth_to_start_2way += od_layer_outward + od_layer_inward;
      }
    }
  }

  /* Ensure that the source power accounts for the two-way optical
     depth to the start of the multiple scattering region */
  config->ss_multiplier *= exp(-optical_depth_to_start_2way);


  if (first_gate < 0 || first_gate >= n-2) {
    /* No gate with small enough mean-free-path found: skip wide-angle
       calculation all together */
    return MS_SUCCESS;
  }

  if (ms_range_spacing_is_regular(n, range, MS_RANGE_SPACING_TOLERANCE)) {
    /* Perform wide-angle calculation on the native grid */
    if (ext_air && ssa_air) {
      return ms_wide_angle(n-first_gate, m-first_gate,
			   config, instrument, surface, range+first_gate,
			   radius+first_gate, ext+first_gate, ssa+first_gate,
			   g+first_gate, ext_air+first_gate, ssa_air+first_gate,
			   bscat_out+first_gate);
    }
    else {
      return ms_wide_angle(n-first_gate, m-first_gate,
			   config, instrument, surface, range+first_gate,
			   radius+first_gate, ext+first_gate, ssa+first_gate,
			   g+first_gate, NULL, NULL,
			   bscat_out+first_gate);
    }
  }
  else {
    /* Regrid the variables before calling ms_wide_angle() */
    ms_real start_edge_range;
    int i;

    /* Calculate the range of the edge of the first multiple
       scattering gate, treating the grid edges as half way between
       the elements of the range vector */ 
    if (first_gate == 0) {
      /* Extrapolate */
      start_edge_range = 1.5*range[0] - 0.5*range[1];
    }
    else {
      /* Interpolate */
      start_edge_range = 0.5*(range[first_gate-1] + range[first_gate]);
    }
    
    /* Find the smallest grid spacing in the first
       MS_OD_DRANGE_THRESHOLD scattering optical depths of the cloud;
       this grid spacing will become "drange_new" */

    /* Far edge of first gate */
    ms_real edge_range2 = 0.5*(range[first_gate] + range[first_gate+1]);

    /* Initial value for drange_new is the spacing for the first gate */
    ms_real drange_new = fabs(start_edge_range - edge_range2);

    /* Initial scattering optical depth */
    ms_real scat_od = drange_new*ext[first_gate]*ssa[first_gate];

    /* Loop through subsequent gates */
    for (i = first_gate+1; i < n; i++) {
      ms_real drange_new_i;     /* Grid spacing of gate i */
      ms_real edge_range1 = edge_range2;/* Near edge of gate i is far edge of
				   gate i-1 */
      /* Calculate location of far edge of gate i */
      if (i < n-1) {
	/* Interpolate */
	edge_range2 = 0.5*(range[i] + range[i+1]);
      }
      else {
	/* Extrapolate */
	edge_range2 = 1.5*range[i] - 0.5*range[i-1];
      }
      /* Set the grid spacing of gate i */
      drange_new_i = fabs(edge_range1-edge_range2);
      /* Is this smaller than drange_new, and is there multiple-scattering
	 cloud present? */
      if (drange_new > drange_new_i && ext[i] > 0.0 && ssa[i] > 0.0) {
	//	ms_real footprint = fabs(instrument.altitude-range[i])
	//	  *instrument.rho_receiver[instrument.nfov-1];
	//	ms_real mfp = 1.0/(ext[i]*ssa[i]*(1.0 - g[i]));
	//	fprintf(stderr, "??? %g %g\n", mfp, footprint);
	//	if (mfp < footprint*MS_MFP_FOOTPRINT_THRESHOLD) {
	  /* Threshold met: set this to the new grid spacing */
	drange_new = drange_new_i;
	  //	}
      }
      /* Increment the scattering optical depth */
      scat_od += drange_new*ext[i]*ssa[i];
      if (scat_od > MS_OD_DRANGE_THRESHOLD) {
	/* We have reached our threshold scattering optical depth: stop
	   looking for more finely spaced gates */
	break;
      }
    }
    /* Calculate location of far edge of the final gate */
    ms_real end_edge_range = 1.5*range[n-1] - 0.5*range[n-2];

    /* Calculate how many points will be required */
    int m = ((int) (0.5+fabs(start_edge_range - end_edge_range)/drange_new));

    /* Allocate intermediate arrays on the stack */
    ms_real range_new[m];
    ms_real radius_new[m];
    ms_real ext_new[m];
    ms_real ssa_new[m];
    ms_real _ext_air_new[m];
    ms_real _ssa_air_new[m];
    ms_real g_new[m];
    ms_real bscat_out_new[m];
    ms_real* ext_air_new = _ext_air_new;
    ms_real* ssa_air_new = _ssa_air_new;
    if (!(ext_air && ssa_air)) {
      ext_air_new = ssa_air_new = NULL;
    }

    int j;
    for (j = 0; j < m; j++) {
      radius_new[j] = 1.0;
      ext_new[j] = ssa_new[j] = g_new[j] = bscat_out_new[j] = 0.0;
      if (ext_air && ssa_air) {
	ext_air_new[j] = ssa_air_new[j] = 0.0;
      }
    }

    /* Interpolate scattering properties on to a regular grid */

    /* Select the first gate of the irregular grid that will be used in
       the interpolation */
    /*
    if (first_gate > 0) {
      i = first_gate-1;
    }
    else {
      i = 0;
    }
    */
    i = first_gate;
    j = 0;
    ms_real direction = (range[1] > range [0] ? 1.0 : -1.0);
    /* Set the range in metres */
    range_new[j] = start_edge_range + direction*drange_new*(0.5+j);
    ms_real edge1_new = range_new[j] - direction*drange_new*0.5;
    ms_real edge2_new = range_new[j] + direction*drange_new*0.5;

    ms_real drange = ms_get_drange(i, n, range);
    //    ms_real edge1 = range[i] - direction*drange*0.5;
    //    ms_real edge2 = range[i] + direction*drange*0.5;
    ms_real edge1 = ms_get_midpoint(i, n, range);
    ms_real edge2 = ms_get_midpoint(i+1, n, range);

    while (j < m) {
      ms_real overlap; /* Length of overlap region in metres */
      if (direction > 0.0) {
	ms_real near = (edge1 > edge1_new ? edge1 : edge1_new);
	ms_real far  = (edge2 > edge2_new ? edge2_new : edge2);
	overlap = far - near;
      }
      else {
	ms_real near = (edge1 < edge1_new ? edge1 : edge1_new);
	ms_real far  = (edge2 < edge2_new ? edge2_new : edge2);
	overlap = near - far;
      }
      ext_new[j] += overlap * ext[i]; /* Optical depth for the
					   moment */
      ssa_new[j] += overlap * ext[i] * ssa[i]; /* Scattering optical
						  depth for the
						  moment */
      g_new[j] += overlap * ext[i] * ssa[i] * g[i];
      radius_new[j] += overlap * ext[i] / (radius[i]*radius[i]);
      if (ext_air && ssa_air) {
	ext_air_new[j] += overlap * ext_air[i];
	ssa_air_new[j] += overlap * ext_air[i] * ssa_air[i];
      }

      ms_real edge2_diff = (edge2_new-edge2)*direction;
      if (edge2_diff <= 0.0 || (edge2_diff >= 0.0 && i >= n-1)) {
	/* Regular grid needs moving forward */
	/* Normalize the values on the regular grid appropriately */
	if (ssa_new[j] > 0.0) {
	  /* At this stage, ssa_new is the scattering optical depth */
	  g_new[j] /= ssa_new[j];
	}
	else {
	  g_new[j] = 0.65;
	}
	if (ext_new[j] > 0.0) {
	  /* At this stage, ext_new is the optical depth */
	  ssa_new[j] /= ext_new[j];
	  /* ...while radius_new is the extinction-weighted integral
	     of radius^-2 */
	  radius_new[j] = sqrt(ext_new[j]/radius_new[j]);
	}
	else {
	  ssa_new[j] = 1.0;
	  radius_new[j] = 1.0;
	}
	ext_new[j] /= drange_new;

	if (ext_air && ssa_air) {
	  if (ext_air_new[j] > 0.0) {
	    ssa_air_new[j] /= ext_air_new[j];
	  }
	  else {
	    ssa_air_new[j] = 0.0;
	  }
	  ext_air_new[j] /= drange_new;
	}

	/* Step the regular grid forward */
	j++;
	if (j >= m) {
	  break;
	}
	range_new[j] = start_edge_range + direction*drange_new*(0.5+j);
	edge1_new = edge2_new;
	edge2_new = range_new[j] + direction*drange_new*0.5;
      }
      if (edge2_diff >= 0.0) {
	i++;
	if (i >= n) {
	  break;
	}
	drange = ms_get_drange(i, n, range);
	edge1 = edge2;
	//	edge2 = range[i] + direction*drange*0.5;
	edge2 = ms_get_midpoint(i+1, n, range);

      }
    }

    /*
    ms_real od = 0.0;
    for (i = 0; i < n; i++) {
      od += ext[i]*ms_get_drange(i, n, range);
    }
    fprintf(stderr, "Old optical depth = %g\n", od);
    
    od = 0.0;
    for (j = 0; j < m; j++) {
      od += ext_new[j]*drange_new;
    }
    fprintf(stderr, "New optical depth = %g\n", od);
    */

    /* Call the TDTS algorithm */
    int status = ms_wide_angle(m, m, config, instrument, surface,
				   range_new, radius_new, ext_new,
				   ssa_new, g_new,
				   ext_air_new, ssa_air_new,
				   bscat_out_new);
    /* Exit if an error occurred */
    if (status != MS_SUCCESS) {
      return status;
    }
    /*
      fprintf(stderr, "  i    range     ext      ssa         g      ext_air      ssa_air     radius    bscat\n");
      for (j = 0; j < m; j++) {
      fprintf(stderr, "%3d %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g\n",
	      j, range_new[j], ext_new[j], ssa_new[j], g_new[j],
	      ext_air_new[j], ssa_air_new[j], radius_new[j], bscat_out_new[j]);
	      }
    */
    
    /* Now interpolate the backscatter back on to the irregular grid,
       but conserving the integral */
    
    i = first_gate;
    j = 0;
    edge1_new = range_new[j] - direction*drange_new*0.5;
    edge2_new = range_new[j] + direction*drange_new*0.5;

    drange = ms_get_drange(i, n, range);
    //    edge1 = range[i] - direction*drange*0.5;
    //    edge2 = range[i] + direction*drange*0.5;
    edge1 = ms_get_midpoint(i, n, range);
    edge2 = ms_get_midpoint(i+1, n, range);

    ms_real bscat = 0.0;
    while (i < n) {
      ms_real overlap; /* Length of overlap region in metres */
      if (direction > 0.0) {
	ms_real near = (edge1 > edge1_new ? edge1 : edge1_new);
	ms_real far  = (edge2 > edge2_new ? edge2_new : edge2);
	overlap = far - near;
      }
      else {
	ms_real near = (edge1 < edge1_new ? edge1 : edge1_new);
	ms_real far  = (edge2 < edge2_new ? edge2_new : edge2);
	overlap = near - far;
      }
      bscat += overlap * bscat_out_new[j];

      ms_real edge2_diff = (edge2_new-edge2)*direction;
      if (edge2_diff >= 0.0 || (edge2_diff <= 0.0 && j >= m-1)) {
	bscat_out[i] += bscat / drange;
	//	fprintf(stderr, "* %g\n", drange);
	bscat = 0.0;
	i++;
	if (i >= n) {
	  break;
	}
	drange = ms_get_drange(i, n, range);
	edge1 = edge2;
	//	edge2 = range[i] + direction*drange*0.5;
	edge2 = ms_get_midpoint(i+1, n, range);

      }
      if (edge2_diff <= 0.0) {
	/* Step the regular grid forward */
	j++;
	if (j >= m) {
	  break;
	}
	edge1_new = edge2_new;
	edge2_new = range_new[j] + direction*drange_new*0.5;
      }
    }
  }
  return MS_SUCCESS;
}

