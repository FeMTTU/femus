/* variance.c -- Fast calculation of the lateral variance of outgoing photons

   Copyright (C) 2006 Robin Hogan <r.j.hogan@reading.ac.uk> 

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

#include <math.h>

#include "ms.h"
#include "ms_switch.h"

/* Calculate the lateral variance of the outgoing unscattered and
   forward-scattered photon distribution */
int
ms_variance(int n, ms_real wavelength, ms_real rho_laser,
	    ms_real inst_altitude, const ms_real *range, 
	    const areal *radius, const areal *ext,
	    areal *variance_out)
{
  areal angle_var = rho_laser*rho_laser;
  areal pos_var = variance_out[0]
    = (range[0]-inst_altitude)*(range[0]-inst_altitude)*angle_var;
  areal covar = fabs(range[0]-inst_altitude)*angle_var;

  int i;
  for (i = 0; i < n-1; i++) {
    ms_real drange = ms_get_drange(i, n, range);
    areal theta = wavelength/(MS_PI*radius[i]);
    pos_var += angle_var*drange*drange + covar*drange;
    covar += angle_var*drange;
    angle_var += 0.5*ext[i]*theta*theta*drange;
    variance_out[i+1] = pos_var;
  }
  return MS_SUCCESS;
}

