/* distribution.c -- Functions to output photon distribution

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

static int msd_num_gates_allocated = 0;
static int msd_num_gates = 0;
static int msd_num_x_allocated = 0;
static int msd_num_x = 0;
static int msd_num_gaussians = 0;
static ms_real msd_dx = 0.0;
static ms_real *msd_Edata;
static ms_real **msd_Edouble;
static ms_real **msd_Emulti;

int
ms_set_distribution_resolution(int num_x, ms_real dx)
{
  if (num_x < 1) {
    return MS_FAILURE;
  }
  else {
    msd_num_x = num_x;
  }
  if (dx <= 0.0) {
    return MS_FAILURE;
  }
  else {
    msd_dx = dx;
  }
  return MS_SUCCESS;
}


int
ms_init_distribution_arrays(int num_gates)
{
  /* Macro for allocating or reallocating memory for a single array */
#define ALLOCATE(target, num_bytes) \
  if ((target = (ms_real*)realloc(target, num_bytes)) == NULL) {	\
      target = NULL; return MS_FAILURE; } 
#define ALLOCATEP(target, num_bytes) \
  if ((target = (ms_real**)realloc(target, num_bytes)) == NULL) {	\
      target = NULL; return MS_FAILURE; } 

  int i, j;
  msd_num_gaussians = 0;
  if (msd_num_x < 1) {
    return MS_FAILURE;
  }
  if (msd_num_gates_allocated < num_gates
      || msd_num_x_allocated < msd_num_x) {
    /* Number of gates allocated is less than number requested (which
       may be zero if this function has not been called): (re)allocate
       memory. */
    msd_num_gates_allocated = 0; /* In case error occurs */
    msd_num_x_allocated = 0;
    ALLOCATE(msd_Edata, 2*num_gates*msd_num_x*sizeof(double));
    ALLOCATEP(msd_Edouble, num_gates*sizeof(double *));
    ALLOCATEP(msd_Emulti, num_gates*sizeof(double *));
    msd_num_gates_allocated = num_gates;
    msd_num_x_allocated = msd_num_x;

    for (i = 0; i < num_gates; i++) {
      msd_Edouble[i] = msd_Edata + i*2*msd_num_x;
      msd_Emulti[i] = msd_Edata + (i*2+1)*msd_num_x;
    }
  }
  for (i = 0; i < num_gates; i++) {
    for (j = 0; j < msd_num_x; j++) {
      msd_Edouble[i][j] = 0.0;
      msd_Emulti[i][j] = 0.0;
    }
  }
  msd_num_gates = num_gates;
  return MS_SUCCESS;
}

int
ms_increment_distribution(int order, int gate, 
			   ms_real energy, ms_real variance)
{
  int i;
  ms_real prefix = energy/(MS_PI*variance);
  ms_real factor = msd_dx*msd_dx/variance;
  /*  fprintf(stderr, "Gate %d, Order %d: %g\n", gate, order, variance); */

  if (order < 2) {
    exit(1);
    return MS_FAILURE;
  }
  if (order == 2) {
    for (i = 0; i < msd_num_x; i++) {
      msd_Edouble[gate][i] += prefix*exp(-factor*(i+0.5)*(i+0.5));
    }
  }
  else {
    for (i = 0; i < msd_num_x; i++) {
      msd_Emulti[gate][i] += prefix*exp(-factor*(i+0.5)*(i+0.5));
    }
  }
  msd_num_gaussians++;
  return MS_SUCCESS;
}

int
ms_print_distribution(ms_config config, FILE* file)
{
  int gate, i;
  for (i = 0; i < msd_num_x; i++) {
    fprintf(file, " %g", msd_dx*(i+0.5));
  }
  fprintf(file, "\n");
  for (gate = 0; gate < msd_num_gates; gate++) {
    for (i = 0; i < msd_num_x; i++) {
      fprintf(file, " %g", msd_Edouble[gate][i]);
    }
    fprintf(file, "\n");
  }
  for (gate = 0; gate < msd_num_gates; gate++) {
    for (i = 0; i < msd_num_x; i++) {
      fprintf(file, " %g", msd_Emulti[gate][i]);
    }
    fprintf(file, "\n");
  }
  if (!(config.options & MS_QUIET)) {
    fprintf(stderr, "Number of Gaussians summed: %d\n", msd_num_gaussians);
  }
  return MS_SUCCESS;
}
