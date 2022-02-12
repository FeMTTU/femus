/* range_spacing_is_regular.c -- Check even spacing of range gates

   Copyright (C) 2004-2007 Robin Hogan <r.j.hogan@reading.ac.uk> 

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
#include <stdio.h>

#include "ms.h"

/* Check that range gate spacing is monotonic and that the spacings
   are the same to within the specified tolerance (e.g. 0.05 for 5% or
   0.01 for 1%; return 1 if regular and 0 if irregular. */
int
ms_range_spacing_is_regular(int n, const ms_real* range, ms_real tolerance)
{
  ms_real max_drange = -1.0e37;
  ms_real min_drange = 1.0e37;
  int i;
  for (i = 0; i < n-1; i++) {
    ms_real drange = range[i+1] - range[i];
    if (drange > max_drange) {
      max_drange = drange;
    }
    if (drange < min_drange) {
      min_drange = drange;
    }
  }

  if (max_drange * min_drange < 0.0) {
    /* Range is not monotonic */
    return 0;
  }

  if (tolerance >= 1.0 || tolerance <= 0.0) {
    fprintf(stderr, "Warning: tolerance argument to ms_check_range_spacing should be greater than 0 and less than 1\n");
  }

  if (fabs(max_drange) > fabs(min_drange)*(1.0+tolerance)
      || fabs(min_drange) > fabs(max_drange)*(1.0+tolerance)) {
    return 0;
  }

  return 1;
}

/* Get the range gate spacing for gate i given n gates with ranges
   centred at "range". The mid-points are assumed to lie half-way
   between the centres, and the end points are treated as one would
   expect. */
ms_real
ms_get_drange(int i, int n, const ms_real *range)
{
  if (i <= 0) {
    return fabs(range[1]-range[0]);
  }
  else if (i >= n-1) {
    return fabs(range[n-1]-range[n-2]);
  }
  else {
    return fabs(range[i+1]-range[i-1])*0.5;
  }
}

/* Get the mid-point between gates i-1 and i. */
ms_real
ms_get_midpoint(int i, int n, const ms_real *range)
{
  if (i <= 0) {
    return 1.5*range[0]-0.5*range[1];
  }
  else if (i >= n) {
    return 1.5*range[n-1]-0.5*range[n-2];
  }
  else {
    return 0.5*(range[i-1]+range[i]);
  }
}

