#include <iostream>
#include <cmath>

/* For details of how to call the multiple scattering algorithms, see
   the comments in multiscatter.h, called from ms.h: */
#include "ms.h"
#include "ms_switch.h"

//  const areal *radius,    /* Cloud/aerosol equivalent radius, microns */
//  const areal *ext,       /* Total extinction coefficient, m-1 */
//  const areal *ssa,       /* Total single-scatter albedo */
//  const areal *g,         /* Total asymmetry factor */
//  const areal *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
//  const ms_real *ext_air,   /* Air ext. coefft., m-1 (NULL for vacuum) */
//  const ms_real *ssa_air,   /* Air single-scatter albedo */
//  const areal *droplet_fraction,/* Fraction of extinction from droplets */
//  const areal *pristine_ice_fraction,/* Fraction of ext from pristine ice */



/* Perform the multiple scattering calculation, using which ever
   combination of algorithms is appropriate depending on multiple
   scattering regime */
int
simple_algorithm(
    /* Input data */
    int n,                   
    int m,                   
    ms_config *config,       
    ms_instrument instrument,
    ms_surface surface,   
    const ms_real *dummy1,
    const areal *a,
    const areal *b,
    const areal *c,
    const areal *d,
    const areal *e,
    const ms_real *dummy2,
    const ms_real *dummy3,
    const areal *f,
    const areal *g,
    /* Output data */
    areal *out,      
    areal *dummy_out)
{

  /*
  for (int j = 0; j < 1000; j++) {
    areal sum = 0.0;
    for (int i = 0; i < n; i++) {
      sum += 2.0*a[i]*exp(b[i])/e[i];
      areal tmp = sum - abs(d[i]-f[i]);
      out[i] = tmp - cos(c[i]*tmp)*(g[i]+1.0);
      //    sum += a[i]*b[i]*c[i];
      //    out[i] = sum + d[i]*e[i]*f[i]*g[i];
    }
  }
  */


  for (int i = 0; i < n; i++) {
    out[i] = 0.0;
  }

  areal product;
  for (int j = 0; j < 1000; j++) {
    product = 1.0;
    for (int i = 0; i < n; i++) {
      product = product * sin(2.0*b[i]+c[i]*d[i]);
    }
    out[0] = out[0] + product;
  }
  

  return MS_SUCCESS;
}
