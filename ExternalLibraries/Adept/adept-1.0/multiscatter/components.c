/* components.c -- Storage of components of small-angle multiscatter calculation

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

/* Single-scattering only backscatter, m-1 sr-1 */
ms_real *ms_sa_bscat_single = NULL;

/* Double-scattering only backscatter, m-1 sr-1 */
ms_real *ms_sa_bscat_double = NULL;

/* Triple and higher scattering, m-1 sr-1 */
ms_real *ms_sa_bscat_multi = NULL;

/* Single-scattering only backscatter, m-1 sr-1 */
ms_real *ms_sa_bscat_air_single = NULL;

/* Double-scattering only backscatter, m-1 sr-1 */
ms_real *ms_sa_bscat_air_double = NULL;

/* Triple and higher scattering, m-1 sr-1 */
ms_real *ms_sa_bscat_air_multi = NULL;

/* Lag of small-angle photons, m */
ms_real *ms_sa_lag;

/* Number of gates currently allocated in the intermediate
   variables: */
static int ms_num_gates_allocated = 0;

/* Number of gates defined in the last call to the algorithm: */
static int ms_num_gates_defined = 0;


/* Intermediate arrays are allocated dynamically and can be read by
   external functions after the algorithm has been run. For speed,
   subsequent runs of the algorithm will use the same memory but new
   memory will be allocated if there is not enough space for the new
   profile. Whether or not new memory is allocated, the arrays are
   always set to zero. Note that this function is called automatically
   by multiscater(). It returns MS_FAILURE if there is a problem
   allocating the memory, MS_SUCCESS otherwise. */
int
ms_init_intermediate_arrays(int num_gates)
{
  /* Macro for allocating or reallocating memory for a single array */
#define ALLOCATE(target, num_bytes) \
   if ((target = realloc(target, num_bytes)) == NULL) { \
      target = NULL; return MS_FAILURE; }

  int i;
  if (ms_num_gates_allocated < num_gates) {
    /* Number of gates allocated is less than number requested (which
       may be zero if this function has not been called): (re)allocate
       memory. */
    int num_bytes = num_gates * sizeof(ms_real);
    ms_num_gates_allocated = 0; /* In case error occurs */
    ALLOCATE(ms_sa_bscat_single, num_bytes);
    ALLOCATE(ms_sa_bscat_double, num_bytes);
    ALLOCATE(ms_sa_bscat_multi, num_bytes);
    ALLOCATE(ms_sa_bscat_air_single, num_bytes);
    ALLOCATE(ms_sa_bscat_air_double, num_bytes);
    ALLOCATE(ms_sa_bscat_air_multi, num_bytes);
    ALLOCATE(ms_sa_lag, num_bytes);
  }
  ms_num_gates_allocated = num_gates;
  /* Set each element of each array to zero. */
  for (i = 0; i < num_gates; i++) {
    ms_sa_bscat_single[i] = 0.0;
    ms_sa_bscat_double[i] = 0.0;
    ms_sa_bscat_multi[i] = 0.0;
    ms_sa_bscat_air_single[i] = 0.0;
    ms_sa_bscat_air_double[i] = 0.0;
    ms_sa_bscat_air_multi[i] = 0.0;
    ms_sa_lag[i] = 0.0;
  }
  ms_num_gates_defined = num_gates;
  return MS_SUCCESS;
}

/* Optional: deallocate memory used by intermediate arrays. */
void
ms_free_intermediate_arrays()
{
  if (ms_num_gates_allocated > 0) {
    free(ms_sa_bscat_single);
    free(ms_sa_bscat_double);
    free(ms_sa_bscat_multi);
    free(ms_sa_bscat_air_single);
    free(ms_sa_bscat_air_double);
    free(ms_sa_bscat_air_multi);
    free(ms_sa_lag);
  }
  ms_num_gates_allocated = 0;
  ms_num_gates_defined = 0;
}


