/* multiscatter_ascii.c -- Single-profile interface for lidar multiple
   scattering algorithm

   Copyright (C) 2004-2011 Robin Hogan <r.j.hogan@reading.ac.uk> 

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


   Run "./multiscatter -help" to get usage information.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#define _GNU_SOURCE 1
#include <fenv.h>

#include "multiscatter.h"

#define MAX_CHARS 128

#ifdef SINGLE_PRECISION
#define READ_HEADER_FORMAT "%d %g %g %g %g\n"
#define READ_FORMAT "%g %g %g %g %g %g %g %g %g %g %n\n"
#else
#define READ_HEADER_FORMAT "%d %lg %lg %lg %lg\n"
#define READ_FORMAT "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %n\n"
#endif

#define CHECK(function) if ((status = (function))) { fprintf(stderr, \
   "Error at line %d of %s: code %d\n", __LINE__, __FILE__, status); exit(status); }

//#define multiscatter simple_algorithm


/* Print usage information */
static
void
usage(char *exec_name)
{
  fprintf(stderr,
	  "Usage\n"
	  "  %s [options] data_in.dat > data_out.dat\n"
	  "  %s [options] < data_in.dat > data_out.dat\n"
	  "\n", exec_name, exec_name);
  fprintf(stderr, 
	  "General options\n" 
	  "  -help           Display this message\n"
	  "  -repeat n       Repeat algorithm n times (for benchmarking)\n"
	  "  -quiet          Don't report activity to stderr\n"
	  "  -v1             Use version 1.x interpretation of first line of input file\n"
	  "  -auto           Automatically select algorithm settings from wavelength etc\n"
	  "  -radar          Use settings appropriate for radar\n"
	  "  -lidar          Use settings appropriate for lidar\n"
	  "  -algorithms <small_angle_algorithm> <wide_angle_algorithm>\n"
	  "                  Manually select algorithm, where <small_angle_algorithm> is:\n"
	  "      none     - No single or small-angle scattering\n"
	  "      single   - Single scattering (no small-angle scattering)\n"
	  "      original - Original Hogan (2006) algorithm: speed O(N^2)\n"
	  "      fast     - Faster Hogan (2008) algorithm: speed O(N) DEFAULT\n"
	  "      explicit - Eloranta-like explicit representation of each order of scattering\n"
	  "      lag      - Hogan (2008) but also with lag calculation (m)\n"
	  "                  ...and where <wide_angle_algorithm> is:\n"
	  "      none     - No wide-angle scattering\n"
	  "      tdts     - Time-dependent two stream (Hogan and Battaglia 2008)\n"
	  "      lidar    - TDTS assuming forward lobe DEFAULT\n"
	  "      radar    - TDTS assuming no forward lobe\n"
	  /* 
	  "  -single-only    Single scattering only\n"
	  "  -sa-only        Don't include wide-angle multiple scattering\n"
	  "  -qsa-only       (alternative to \"-sa-only\")\n"
	  "  -wide-only      Only wide-angle multiple scattering\n"
	  */
	  "  -hsrl           Output particulate and air backscatter separately\n"
	  "  -gaussian-receiver\n"
	  "                  Receiver is Gaussian rather than top-hat shaped\n"
	  "  -jacobian       Output the approximate but fast Jacobian\n"
	  "  -numerical-jacobian\n"
	  "                  Output the Jacobian calculated (slowly) using finite\n"
	  "                  differences\n"
	  "  -ext-only       Only calculate the Jacobian with respect to extinction\n"
	  "  -adjoint        Output the adjoint as well\n"
	  "  -check-adjoint  Calculate the adjoint and check it with a numerical Jacobian\n"
	  "  -annular-detectors\n"
	  "                  2nd and subsequent detectors are ring shaped\n"
	  "\n");
  fprintf(stderr,
	  "Options for ORIGINAL small-angle (SA) algorithm\n"
	  /*
	  "  -fast-sa        Use fast O(N) SA model\n"
	  "  -fast-qsa       (alternative to \"-fast-sa\")\n"
	  "  -lag            Calculate SA lag (m) - EXPERIMENTAL\n"
	  */
	  "  -simple-optical-depth\n"
	  "                  Use simple optical depth integration\n"
	  "  -crude-optical-depth\n"
	  "                  Use crude optical depth integration\n"
	  "  -crude-integration\n"
	  "                  Use crude double/multiple scattering integration\n");
  fprintf(stderr,
	  "  -no-multiscat-within-gate\n"
	  "                  Photon cannot be forward scattered more than once within a\n"
	  "                  single gate; note that this has the opposite effect to\n"
	  "                  \"-crude-integration\"\n"
	  "  -double-scattering-only\n"
	  "                  Don't include triple and higher-order scatterings\n");
  fprintf(stderr,
	  "  -no-molecular-extinction\n"
	  "                  Molecules do not extinguish the beam\n"
	  "  -wide-angle-cutoff <theta>\n"
	  "                  Forward scattering at angles greater than <theta> radians\n"
	  "                  are deemed to escape, a crude way to deal with a problem\n"
	  "                  associated with aerosols\n");
  fprintf(stderr,
	  "  -crude-double-scattering\n"
	  "                  Don't use Eloranta's slow but accurate double scattering\n"
	  "                 formulation\n"
	  "  -approx-exp     Appriximate the exp() function call for speed\n"
	  "\n");
  fprintf(stderr,
	  "Options for EXPLICIT small-angle (SA) algorithm\n"
	  "  -explicit n     Use an explicit model with n orders of scattering\n"
	  "  -output-distribution n dx\n"
	  "                  Output horizontal photon distributions at n points spaced\n"
	  "                  dx apart, starting at dx/2\n"
	  "\n");
  fprintf(stderr,
	  "Options for wide-angle multiple-scattering\n"
	  /*
	  "  -no-forward-lobe\n"
	  "                  Radar-like phase function behaviour: use single-scattering\n"
	  "                  rather than SA\n"
	  */
	  "  -ignore-source-widening\n"
	  "                  Ignore widening effect of small-angle scattering on source beam for wide-angle calc\n"
	  "  -optimize-wide-angle-gates\n"
	  "                  Skip optically thin gates in wide-angle calculation\n"
	  "  -simple-2s-coeffts\n"
	  "                  Use the simple upwind Euler formulation (can be unstable\n"
	  "                  for high optical depth)\n");
  fprintf(stderr,
	  "  -ssa-scales-forward-lobe\n"
	  "                  Single-scattering albedo less than unity reduces the\n"
	  "                  forward scattering lobe as well as wide-angle scattering\n"
	  "  -num-samples m  Output m samples, allowing sampling of signals appearing\n"
	  "                  to originate below ground\n"
	  "  -propagation-to-stderr\n"
	  "                  Print the diffuse photon energies and variances to\n"
	  "                  standard error\n"
	  "\n");
  fprintf(stderr,
	  "Input data\n"
	  "  First line: 5 values\n"
	  "    1: Number of range gates in profile\n"
	  "    2: Wavelength (m)\n"
	  "    3: Altitude of instrument (m)\n"
	  "    4: Transmitter divergence, 1/e half-width (radians)\n"
	  "    5: Receiver field-of-view, half-width of first receiver (radians)\n"
	  "   6+: Fields-of-view for any other channels (radians)\n");
  fprintf(stderr,
	  "  (note that with the \"-v1\" option, the version 1.x file format is assumed,\n"
	  "   in which the order is (1) number of range gates, (2) wavelength, (3) divergence,\n"
	  "   (4) field-of-view, (5) altitude, and only one field of view is allowed)\n");
  fprintf(stderr,
	  "  Subsequent lines: 4, 5, 8, 10 or more values (all 8 required for the wide-angle\n"
	  "  calculation; all 10 to represent anisotropic phase functions near 180 deg);\n"
	  "  more to calculate the adjoint\n"
	  "    1: Range of gate above ground starting with nearest gate to instrument (m)\n"
	  "    2: Extinction coefficient of cloud/aerosol only (m-1)\n"
	  "    3: Equivalent-area cloud/aerosol radius (m)\n"
	  "    4: Extinction-to-backscatter ratio of cloud/aerosol (sr)\n"
	  "    5: Air extinction coefficient (m-1); optional (zero if missing)\n");
  fprintf(stderr,
	  "    6: Single scattering albedo of cloud/aerosol\n"
	  "    7: Scattering asymmetry factor of cloud/aerosol\n"
	  "    8: Single scattering albedo of air\n");
  fprintf(stderr,
	  "    9: Fraction of cloud/aerosol backscatter due to droplets (real number\n"
	  "         between 0 and 1), for representing anisotropic near-180 phase\n"
	  "         functions in the SA algorithm\n"
	  "   10: Pristine ice fraction, for anisotropic near-180 phase functions\n"
	  "         from Yang et al. (2000); note that irregular ice particles\n"
	  "         are believed to have isotropic near-180 phase functions, in\n"
	  "         which case set this variable to zero\n");
  fprintf(stderr, 
	  "  Subsequent lines (only with \"-adjoint\" option)\n"
	  " 11 to 10+n: Adjoint input for particulate backscatter of the n fields-of-view\n"
	  " 11+n: Adjoint input for air backscatter of the first field-of-view (for HSRL)\n"
	  "\n");
  fprintf(stderr,
	  "Output data (without -jacobian or -numerical-jacobian option)\n"
	  "  Default: 5 values per range gate\n"
	  "    1: Index of range gate\n"
	  "    2: Apparent range above ground (m)\n"
	  "    3: Extinction coefficient (m-1)\n"
	  "    4: Equivalent-area particle radius (microns)\n"
	  "    5: Apparent backscatter (m-1 sr-1), for particulates only if with \"-hsrl\"\n");
  fprintf(stderr,
	  "  With the \"-hsrl\" option\n"
	  "    6: Apparent backscatter of the air only (m-1 sr-1)\n");
  fprintf(stderr,
	  "  With the \"-adjoint\" option (where n=6 if \"-hsrl\", n=5 otherwise)\n"
	  "  n+1: Adjoint of extinction coefficient (m)\n"
	  "  n+2: Adjoint of single-scattering albedo\n"
	  "  n+3: Adjoint of asymmetry factor\n"
	  "  n+4: Adjoint of extinction-to-backscatter ratio (m sr)\n"
	  "\n");
  fprintf(stderr,
	  "Output data (with -jacobian or -numerical-jacobian option)\n"
	  "   n x m matrix: d(attenuated backscatter) / d(extinction coefficient)\n"
	  "and if \"-ext-only\" is not set:\n"
	  "   n x m matrix: d(attenuated backscatter) / d(single-scattering albedo)\n"
	  "   n x m matrix: d(attenuated backscatter) / d(asymmetry factor)\n"
	  "   n x m matrix: d(attenuated backscatter) / d(particle radius)\n"
	  "   1 x m vector: d(attenuated backscatter) / d(ext-to-bscat ratio)\n");
}

static
void
print_error_AD(FILE* file, ms_real AD, ms_real AD_test) {
#ifndef PRINT_RAW_AD_ERROR
  if (AD == 0.0) {
    if (AD_test == 0.0) {
      fprintf(file, "   0      ");
    }
    else {
      fprintf(file, "%+g\t", AD_test);
    }
  }
  else {
    fprintf(file, "%+8.3f%% ", 100.0*(AD_test-AD)/AD);
  }
#else
  fprintf(file, "(%g,%g)\t", AD, AD_test);

#endif
}


/* Main program */
int
main(int argc, char **argv)
{
  /* VARIABLE DECLARATIONS */

  FILE *infile = stdin;
  int n, m = 0, i, iarg;
  int nrepeats = 1;
  int output_jacobian = 0;
  int output_adjoint = 0;
  int calc_jacobian = 0;
  int separate_bscat_air = 0;
  int use_isotropic_pp = 0;
  int use_air_ext = 1;
  int automatically_configure = 0;
  int automatically_configure_tdts = 0;
  int manually_select_algorithms = 0;
  int manually_select_cbh = 0;
  int print_stats = 0;
  int use_version_1x = 0;
  int ninputs;
  int norder = 4;
  char buffer[MAX_CHARS];
  int ch, status;
  char* strbuf = NULL;
  char* curbuf = NULL;
  char* newcurbuf = NULL;
  size_t strbufoffset = 0;
  size_t strbuflen = 0;
  ms_real rho_receiver_old = 1.0;

  ms_real* range = NULL;
  ms_real* radius = NULL;
  ms_real* ext = NULL;
  ms_real* ext_bscat_ratio = NULL;
  ms_real* ext_air = NULL;
  ms_real* droplet_fraction = NULL;
  ms_real* pristine_ice_fraction = NULL;
  ms_real* ssa = NULL;
  ms_real* g = NULL;
  ms_real* ssa_air = NULL;

  ms_real* bscat = NULL;
  ms_real* bscat_air = NULL;

  ms_real* ext_AD = NULL;
  ms_real* ssa_AD = NULL;
  ms_real* g_AD = NULL;
  ms_real* radius_AD = NULL;
  ms_real* ext_bscat_ratio_AD = NULL;
  
  ms_real* bscat_AD = NULL;
  ms_real* bscat_air_AD = NULL;

  ms_config config = MS_DEFAULT_CONFIG;
  ms_instrument instrument = MS_DEFAULT_INSTRUMENT;
  ms_surface surface = MS_DEFAULT_SURFACE;

  /* Enable some exceptions. At startup all exceptions are masked. */
   //   feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW|FE_UNDERFLOW);
   //FE_DIVBYZERO, FE_INEXACT, FE_INVALID, FE_OVERFLOW, FE_UNDERFLOW 
   //  feenableexcept(FE_ALL_EXCEPT);
  //  fesetround(FE_UPWARD);
  /* HANDLE COMMAND-LINE ARGUMENTS */

  for (iarg = 1; iarg < argc; ++iarg) {
    if (strcmp(argv[iarg], "-auto") == 0) {
      automatically_configure = 1;
    }
    else if (strcmp(argv[iarg], "-algorithms") == 0) {
      if (iarg + 2 < argc) {
	iarg++;
	manually_select_algorithms = 1;
	if (strcmp(argv[iarg], "none") == 0) {
	  config.small_angle_algorithm = MS_SINGLE_AND_SMALL_ANGLE_NONE;
	}
	else if (strcmp(argv[iarg], "single") == 0) {
	  config.small_angle_algorithm = MS_SINGLE_SCATTERING;
	}
	else if (strcmp(argv[iarg], "original") == 0) {
	  config.small_angle_algorithm = MS_SMALL_ANGLE_PVC_ORIGINAL;
	}
	else if (strcmp(argv[iarg], "fast") == 0) {
	  config.small_angle_algorithm = MS_SMALL_ANGLE_PVC_FAST;
	}
	else if (strcmp(argv[iarg], "explicit") == 0) {
	  config.small_angle_algorithm = MS_SMALL_ANGLE_PVC_EXPLICIT;
	}
	else if (strcmp(argv[iarg], "lag") == 0) {
	  config.small_angle_algorithm = MS_SMALL_ANGLE_PVC_FAST_LAG;
	}
	else {
	  fprintf(stderr, "Error: small angle algorithm \"%s\" not recognised\n",
		  argv[iarg]);
	  exit(MS_COMMAND_LINE_ERROR);
	}
	iarg++;
	if (strcmp(argv[iarg], "none") == 0) {
	  config.wide_angle_algorithm = MS_WIDE_ANGLE_NONE;
	}
	else if (strcmp(argv[iarg], "tdts") == 0) {
	  /* Use the wavelength to decide which version of TDTS is the
	     most appropriate to ensure no energy is lost */
	  automatically_configure_tdts = 1;
	  /*
	  if (config.small_angle_algorithm <= MS_SINGLE_SCATTERING) {
	    config.wide_angle_algorithm = MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE;
	  }
	  else {
	    config.wide_angle_algorithm = MS_WIDE_ANGLE_TDTS_FORWARD_LOBE;
	  }
	  */
	}
	else if (strcmp(argv[iarg], "lidar") == 0) {
	  config.wide_angle_algorithm = MS_WIDE_ANGLE_TDTS_FORWARD_LOBE;
	}
	else if (strcmp(argv[iarg], "radar") == 0) {
	  config.wide_angle_algorithm = MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE;
	}
	else {
	  fprintf(stderr, "Error: wide angle algorithm \"%s\" not recognised\n",
		  argv[iarg]);
	  exit(MS_COMMAND_LINE_ERROR);
	}
      }
      else {
	fprintf(stderr, "Error: option \"-algorithms\" requires the small- and wide-angle algorithms to be listed\n");
	exit(MS_COMMAND_LINE_ERROR);
      }
    }
    else if (strcmp(argv[iarg], "-lidar") == 0) {
      config.small_angle_algorithm = MS_SMALL_ANGLE_PVC_FAST;
      config.wide_angle_algorithm = MS_WIDE_ANGLE_TDTS_FORWARD_LOBE;
    }
    else if (strcmp(argv[iarg], "-radar") == 0) {
      config.small_angle_algorithm = MS_SINGLE_SCATTERING;
      config.wide_angle_algorithm = MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE;
    }
    else if (strcmp(argv[iarg], "-v1") == 0) {
      /* Use version 1.x file format */
      use_version_1x = 1;
    }
    else if (strcmp(argv[iarg], "-hsrl") == 0) {
      /* Output separate particulate and air backscatters */
      separate_bscat_air = 1;
    }
    else if (strcmp(argv[iarg], "-annular-detectors") == 0) {
      /* Output separate particulate and air backscatters */
      config.options |= MS_ANNULAR_DETECTORS;
    }
    //    else if (strcmp(argv[iarg], "-lag") == 0) {
      /* Use fast O(N) small-angle calculation with lag
	 calculation */
    //      config.options |= MS_FAST_QSA;
    //      config.options |= MS_QSA_LAG;
    //    }
    else if (strcmp(argv[iarg], "-output-bscats") == 0) {
      /* Output different backscatters */
      fprintf(stderr, "Error: the %s option is deprecated\n", argv[iarg]);
      exit(MS_COMMAND_LINE_ERROR);
    }
    else if (strcmp(argv[iarg], "-output-all") == 0) {
      /* Output all internal variables */
      fprintf(stderr, "Error: the %s option is deprecated\n", argv[iarg]);
      exit(MS_COMMAND_LINE_ERROR);
    }
    else if (strcmp(argv[iarg], "-output-distribution") ==0) {
      int n;
      ms_real dx;
      if (argc > ++iarg) {
	n = atoi(argv[iarg]);
	if (n < 1) {
	  fprintf(stderr,
		  "Error: number of distribution points must be greater than 0\n");
	  exit(MS_COMMAND_LINE_ERROR);
	}
	if (argc > ++iarg) {
	  dx = atof(argv[iarg]);
	  if (dx < 0.0) {
	    fprintf(stderr,
		    "Error: distribution point spacing must be greater than 0.0\n");
	    exit(MS_COMMAND_LINE_ERROR);
	  }
	  else {
	    config.options |= MS_OUTPUT_DISTRIBUTION;
	    ms_set_distribution_resolution(n, dx);
	  }
	}
	else {
	  fprintf(stderr, "Error: distribution point spacing (dx) not specified\n");
	  exit(MS_COMMAND_LINE_ERROR);
	}
      }
      else {
	fprintf(stderr, "Error: number of distribution points not specified\n");
	exit(MS_COMMAND_LINE_ERROR);
      }
    }
    else if (strcmp(argv[iarg], "-print-stats") == 0) {
      /* Print statistics */
      print_stats = 1;
    }
    else if (strcmp(argv[iarg], "-gaussian-receiver") == 0) {
      /* Use radar-like Gaussian receiver pattern rather than
	 top-hat */
      instrument.receiver_type = MS_GAUSSIAN;
    }
    else if (strcmp(argv[iarg], "-coherent-enhancement") == 0) {
      if (argc > ++iarg) {
	manually_select_cbh = 1;
	float val = atof(argv[iarg]);
	if (val < 1.0 || val > 2.0) {
	  fprintf(stderr,
		  "Error: coherent backscatter enhancement must be between 1.0 and 2.0\n");
	  exit(MS_COMMAND_LINE_ERROR);
	}
	config.coherent_backscatter_enhancement = val;
      }
      else {
	fprintf(stderr, "Error: coherent backscatter enhancement not specified\n");
	exit(MS_COMMAND_LINE_ERROR);
      }      
    }
    else if (strcmp(argv[iarg], "-no-forward-lobe") == 0) {
      /* Use radar-like phase function behaviour: no narrow forward
	 lobe */
      config.options |= MS_NO_FORWARD_LOBE;
    }
    else if (strcmp(argv[iarg], "-ssa-scales-forward-lobe") == 0) {
      config.options |= MS_SSA_SCALES_FORWARD_LOBE;
    }
    else if (strcmp(argv[iarg], "-simple-2s-coeffts") == 0) {
      config.options |= MS_SIMPLE_2S_COEFFTS;
    }
    else if (strcmp(argv[iarg], "-simple-optical-depth") == 0) {
      /* Use simple optical depth calculation */
      config.options |= MS_SIMPLE_OPTICAL_DEPTH;
    }
    else if (strcmp(argv[iarg], "-crude-optical-depth") == 0) {
      /* Use crude optical depth calculation */
      config.options |= MS_CRUDE_OPTICAL_DEPTH;
    }
    else if (strcmp(argv[iarg], "-crude-integration") == 0) {
      /* Use crude integration for double and multiple scattering */
      config.options |= MS_CRUDE_INTEGRATION;
      config.options &= ~MS_NO_MULTISCAT_WITHIN_GATE;
    }
    else if (strcmp(argv[iarg], "-no-multiscat-within-gate") == 0) {
      /* Each photon cannot be forward scattered more than once within
	 a range gate */
      config.options |= MS_NO_MULTISCAT_WITHIN_GATE;
      config.options &= ~MS_CRUDE_INTEGRATION;
    }
    else if (strcmp(argv[iarg], "-no-molecular-extinction") == 0) {
      /* Molecules backscatter but don't extinguish */
      config.options |= MS_NO_MOLECULAR_EXTINCTION;
    }
    else if (strcmp(argv[iarg], "-crude-double-scattering") == 0) {
      /* Use crude but fast method for double scattering */
      config.options |= MS_CRUDE_DOUBLE_SCATTERING;
    }
    else if (strcmp(argv[iarg], "-approx-exp") == 0) {
      /* Use crude but very fast method for double scattering */
      config.options |= MS_APPROXIMATE_EXPONENTIAL;
    }
    else if (strcmp(argv[iarg], "-double-scattering-only") == 0) {
      /* Use crude but fast method for double scattering */
      config.options |= MS_DOUBLE_SCATTERING_ONLY;
    }
    else if (strcmp(argv[iarg], "-ignore-source-widening") == 0) {
      /* Ignore the widening effect of foward scattering on the
	 quasi-direct outgoing beam for the purposes of wide-angle
	 scattering */
      config.options |= MS_IGNORE_SOURCE_WIDENING;
    }
    else if (strcmp(argv[iarg], "-optimize-wide-angle-gates") == 0) {
      /* Ignore the widening effect of foward scattering on the
	 quasi-direct outgoing beam for the purposes of wide-angle
	 scattering */
      config.first_wide_angle_gate = MS_AUTO_FIRST_WIDE_ANGLE_GATE;
    }
    else if (strcmp(argv[iarg], "-wide-angle-cutoff") == 0) {
      /* A crude method to cope with the problem of representing both
	 cloud and wide-angle-scattering aerosol in the same profile,
	 without running the full wide-angle scattering code */
      if (iarg + 1 < argc) {
	config.options |= MS_WIDE_ANGLE_CUTOFF;
	config.max_theta = atof(argv[++iarg]);
      }
      else {
	fprintf(stderr, "Error: option \"-wide-angle-cutoff\" requires an angle to be specified\n");
	exit(MS_COMMAND_LINE_ERROR);
      }
    }
    else if (strcmp(argv[iarg], "-explicit") == 0) {
      if (iarg + 1 < argc) {
	norder = atoi(argv[++iarg]);
	config.small_angle_algorithm = MS_SMALL_ANGLE_PVC_EXPLICIT;
	if (ms_set_max_scattering_order(&config, norder) == MS_FAILURE) {
	  fprintf(stderr, "Error: number of scattering orders in explicit calculation out of range\n");
	  exit(MS_COMMAND_LINE_ERROR);
	}
      }
      else {
	fprintf(stderr, "Error: option \"-explicit\" requires the number of scattering orders to be specified\n");
	exit(MS_COMMAND_LINE_ERROR);
      }
    }
    else if (strcmp(argv[iarg], "-jacobian") == 0) {
      /* Output the Jacobian */
      output_jacobian = 1;
      output_adjoint = 0;
      if (!calc_jacobian) {
	calc_jacobian = 1;
      }
    }
    else if (strcmp(argv[iarg], "-analytic-jacobian") == 0) {
      config.options |= MS_JACOBIAN;
    }
    else if (strcmp(argv[iarg], "-numerical-jacobian") == 0) {
      /* Output the Jacobian */
      output_jacobian = 1;
      output_adjoint = 0;
      if (!calc_jacobian) {
	calc_jacobian = 1;
      }
      config.options |= MS_NUMERICAL_JACOBIAN;
    }
    else if (strcmp(argv[iarg], "-ext-only") == 0) {
      /* Output the Jacobian with respect to extinction  */
      calc_jacobian = 2;
      config.options |= MS_NUMERICAL_JACOBIAN;
    }
    else if (strcmp(argv[iarg], "-adjoint") == 0) {
      output_adjoint = 1;
      output_jacobian = 0;
    }
    else if (strcmp(argv[iarg], "-check-adjoint") == 0) {
      output_adjoint = 2;
      if (!calc_jacobian) {
	calc_jacobian = 1;
      }
      output_jacobian = 0;
    }
    else if (strcmp(argv[iarg], "-adept") == 0) {
      config.options |= MS_AUTOMATIC_DIFFERENTIATION_ADEPT;
    }
    else if (strcmp(argv[iarg], "-cppad") == 0) {
      config.options |= MS_AUTOMATIC_DIFFERENTIATION_CPPAD;
    }
    else if (strcmp(argv[iarg], "-adolc") == 0) {
      config.options |= MS_AUTOMATIC_DIFFERENTIATION_ADOLC;
    }
    else if (strcmp(argv[iarg], "-sacado") == 0) {
      config.options |= MS_AUTOMATIC_DIFFERENTIATION_SACADO;
    }
    else if (strcmp(argv[iarg], "-sacado-fad") == 0) {
      config.options |= MS_AUTOMATIC_DIFFERENTIATION_SACADO_FAD;
    }
    else if (strcmp(argv[iarg], "-force-forward-jacobian") == 0) {
      config.force_jacobian = 1;
    }
    else if (strcmp(argv[iarg], "-force-reverse-jacobian") == 0) {
      config.force_jacobian = -1;
    }
    else if (strcmp(argv[iarg], "-propagation-to-stderr") == 0) {
      config.options |= MS_PROPAGATION_TO_STDERR;
      config.options |= MS_QUIET;
    }
    else if (strcmp(argv[iarg], "-repeat") == 0) {
      /* Perform repeat calculations for benchmarking purposes */
      if (argc > ++iarg) {
	nrepeats = atoi(argv[iarg]);
	if (nrepeats < 1 || nrepeats > 1000000) {
	  fprintf(stderr,
		  "Error: number of repeats must be between 1 and 1000000\n");
	  exit(MS_COMMAND_LINE_ERROR);
	}
      }
      else {
	fprintf(stderr, "Error: number of repeats not specified\n");
	exit(MS_COMMAND_LINE_ERROR);
      }
    }
    else if (strcmp(argv[iarg], "-num-samples") == 0) {
      if (argc > ++iarg) {
	m = atoi(argv[iarg]);
	if (m < 1) {
	  fprintf(stderr,
		  "Error: number of instrument samples must be greater than 0\n");
	  exit(MS_COMMAND_LINE_ERROR);
	}
      }
      else {
	fprintf(stderr, "Error: number of instrument samples not specified\n");
	exit(MS_COMMAND_LINE_ERROR);
      }      
    }
    else if (strcmp(argv[iarg], "-quiet") == 0) {
      config.options |= MS_QUIET;
    }
    else if (strcmp(argv[iarg], "-help") == 0) {
      usage(argv[0]);
      exit(MS_COMMAND_LINE_ERROR);
    }
    else if (strcmp(argv[iarg], "-fast-sa") == 0
	     || strcmp(argv[iarg], "-fast-qsa") == 0
	     || strcmp(argv[iarg], "-lag") == 0
	     || strcmp(argv[iarg], "-no-forward-lobe") == 0
	     || strcmp(argv[iarg], "-single-only") == 0
	     || strcmp(argv[iarg], "-sa-only") == 0
	     || strcmp(argv[iarg], "-qsa-only") == 0
	     || strcmp(argv[iarg], "-wide-only") == 0) {
      fprintf(stderr, "Error: \"%s\" no longer used to control algorithms: use \"-algorithms\" instead\n",
	      argv[iarg]);
      usage(argv[0]);
      exit(MS_COMMAND_LINE_ERROR);
    }
    else if (argv[iarg][0] == '-') {
      /* Argument not understood */
      fprintf(stderr, "Error: \"%s\" not understood\n", argv[iarg]);
      usage(argv[0]);
      exit(MS_COMMAND_LINE_ERROR);
    }
    else {
      /* Assume the argument is a filename */
      if (! (infile = fopen(argv[iarg], "r"))) {
	fprintf(stderr, "Error: \"%s\" not found\n", argv[iarg]);
	exit(MS_COMMAND_LINE_ERROR);
      }
      break;
    }
  }


  /* REPORT WHAT IS BEING DONE */

  if (!(config.options & MS_QUIET)) {
    fprintf(stderr, "Multiscatter %s\n", MS_VERSION);
    fprintf(stderr, "Type \"%s -help\" for usage information\n", argv[0]);
    if (infile == stdin) {
      fprintf(stderr, "Reading stdin...\n");
    }
    else {
      fprintf(stderr, "Reading %s...\n", argv[iarg]);
    }
  }

  /* READ AND CHECK THE INPUT FILE */

  /* Skip commented lines */
  while ((ch = fgetc(infile)) == '#') {
    while ((ch = fgetc(infile)) != '\n') {
      if (ch == EOF) {
	fprintf(stderr, "Error: only comments found in input file\n");
	exit(MS_INPUT_FILE_ERROR);
      }
    }
  }
  ungetc(ch, infile);

  /* Read input data */
  /* First line can be any length, in principle: read it into strbuf */
  do {
    strbuflen += MAX_CHARS;
    strbuf = (char*)realloc(strbuf, sizeof(char)*strbuflen);
    if (!strbuf) {
      fprintf(stderr, "Error allocating memory to read first line of data\n");
      exit(MS_MEMORY_ALLOCATION_ERROR);
    }
    if (!fgets(strbuf+strbufoffset, MAX_CHARS, infile)) {
      fprintf(stderr, "Error reading first line of data\n");
      exit(MS_INPUT_FILE_ERROR);
    }
    /* Search for a new line */
    while (strbuf[strbufoffset] != '\0'
	   && strbuf[strbufoffset] != '\n') {
	/* Found end of line */
      ++strbufoffset;
    }
  }
  while (strbuf[strbufoffset] != '\n');

  /* Now read it in terms of header variables */
  if (use_version_1x) {
    /* Version 1.x file format */
    if (!(config.options & MS_QUIET)) {
      fprintf(stderr, "Assuming version 1.x file header format\n");
    }

    ninputs = fscanf(infile, READ_HEADER_FORMAT, &n,
		     &instrument.wavelength,
		     &instrument.rho_transmitter,
		     &rho_receiver_old, &instrument.altitude);
    if (ninputs < 5) {
      fprintf(stderr, "Error: garbled input at first line: "
	      "only %d inputs read successfully\n", ninputs);
      exit(MS_INPUT_FILE_ERROR);
    }
    instrument.rho_receiver = &rho_receiver_old;
    instrument.nfov = 1;
  }
  else {
    /* Version 2.x file format */
    n = strtol(strbuf, &curbuf, 10);
    instrument.wavelength = strtod(curbuf, &newcurbuf);
    if (instrument.wavelength <= 0.0 || curbuf == newcurbuf) {
      fprintf(stderr, "Error reading wavelength as a positive real\n");
      exit(MS_INPUT_FILE_ERROR);
    }
    
    curbuf = newcurbuf;
    instrument.altitude = strtod(curbuf, &newcurbuf);
    if (curbuf == newcurbuf) {
      fprintf(stderr, "Error reading instrument altitude\n");
      exit(MS_INPUT_FILE_ERROR);
    }
    
    curbuf = newcurbuf;
    instrument.rho_transmitter = strtod(curbuf, &newcurbuf);
    if (instrument.rho_transmitter < 0.0 || curbuf == newcurbuf) {
      fprintf(stderr, "Error reading beam divergence as non-negative real\n");
      exit(MS_INPUT_FILE_ERROR);
    }
    
    /*
      curbuf = newcurbuf;
      instrument.rho_receiver = realloc(instrument.rho_receiver,
      sizeof(ms_real));
      if (!instrument.rho_receiver) {
      fprintf(stderr, "Error allocating data for first receiver field-of-view\n");
      exit(MS_MEMORY_ALLOCATION_ERROR);
      }
      instrument.rho_receiver[0] = strtod(curbuf, &newcurbuf);
      if (curbuf == newcurbuf) {
      fprintf(stderr, "Error reading first field-of-view\n");
      exit(MS_INPUT_FILE_ERROR);
      }
    */
    
    instrument.nfov = -1; /* OK, it should be 1 by now but will soon be
			     incremented... */
    while (curbuf != newcurbuf) {
      curbuf = newcurbuf;
      instrument.nfov++;
      instrument.rho_receiver = (ms_real*)realloc(instrument.rho_receiver,
					sizeof(ms_real)*(instrument.nfov+1));
      if (!instrument.rho_receiver) {
	fprintf(stderr, "Error allocating data for receiver fields-of-view\n");
	exit(MS_MEMORY_ALLOCATION_ERROR);
      }
      instrument.rho_receiver[instrument.nfov] = strtod(curbuf, &newcurbuf);
      if (curbuf != newcurbuf
	  && instrument.rho_receiver[instrument.nfov] <= 0.0) {
	fprintf(stderr, "Error: receiver half-angle field-of-view is %g radians, must be > 0\n",
		instrument.rho_receiver[instrument.nfov]);
	exit(MS_INPUT_FILE_ERROR);
      }
    }
    if (instrument.nfov <= 0) {
      fprintf(stderr, "Error reading receiver fields-of-view\n");
      exit(MS_INPUT_FILE_ERROR);
    }
    free(strbuf);
    
    if (instrument.rho_receiver[0] >= 1.0) {
      fprintf(stderr, "*** WARNING ***\n"
	      "  First receiver field-of-view is %g radians.\n"
	      "  Are you sure you are not using the input format from multiscatter 1.x\n"
	      "  and this is altitude?  If so, either change the top line of the input\n"
	      "  file to use the new format (run \"%s -help\" to see that) or use the\n"
	      "  \"-v1\" option to use the old input format.\n", instrument.rho_receiver[0],
	      argv[0]);
    }
  }

  /* Allocate memory */
  range = (ms_real*)malloc(sizeof(ms_real)*n);
  radius = (ms_real*)malloc(sizeof(ms_real)*n);
  ext = (ms_real*)malloc(sizeof(ms_real)*n);
  ext_bscat_ratio = (ms_real*)malloc(sizeof(ms_real)*n);
  ext_air = (ms_real*)malloc(sizeof(ms_real)*n);
  ssa = (ms_real*)malloc(sizeof(ms_real)*n);
  ssa_air = (ms_real*)malloc(sizeof(ms_real)*n);
  droplet_fraction = (ms_real*)malloc(sizeof(ms_real)*n);
  pristine_ice_fraction = (ms_real*)malloc(sizeof(ms_real)*n);
  g = (ms_real*)malloc(sizeof(ms_real)*n);

  if (!range || !radius || !ext || !ext_bscat_ratio || !ext_air 
      || !ssa || !ssa_air || !droplet_fraction || !pristine_ice_fraction || !g) {
    fprintf(stderr, "Error allocating space for input variables\n");
    exit(MS_MEMORY_ALLOCATION_ERROR);
  }


  /* If number of samples not specified then it is set equal to the
     number of data levels in the profile */
  if (m == 0) {
    m = n;
  }


  /* Allocate space for adjoints if required */
  if (output_adjoint) {

    bscat_AD = (ms_real*)malloc(sizeof(ms_real)*instrument.nfov*m);
    if (separate_bscat_air) {
      bscat_air_AD = (ms_real*)malloc(sizeof(ms_real)*instrument.nfov*m);
    }
    if (!bscat_AD || ((!bscat_air_AD) && separate_bscat_air)) {
      fprintf(stderr, "Error allocating space for adjoint inputs\n");
      exit(MS_MEMORY_ALLOCATION_ERROR);
    }

    ext_AD = (ms_real*)malloc(sizeof(ms_real)*n);
    ssa_AD = (ms_real*)malloc(sizeof(ms_real)*n);
    g_AD = (ms_real*)malloc(sizeof(ms_real)*n);
    radius_AD = (ms_real*)malloc(sizeof(ms_real)*n);
    ext_bscat_ratio_AD = (ms_real*)malloc(sizeof(ms_real)*n);
    if (!ext_AD || !ssa_AD || !g_AD || !ext_bscat_ratio_AD
	|| !radius_AD) {
      fprintf(stderr, "Error allocating space for adjoint outputs\n");
      exit(MS_MEMORY_ALLOCATION_ERROR);
    }
    for (i = 0; i < n; i++) {
      ext_AD[i] = ssa_AD[i] = g_AD[i] = radius_AD[i]
	= ext_bscat_ratio_AD[i] = 0.0;
    }
  }


  /* Read in the n data levels */
  for (i = 0; i < n; i++) {
    int nchars = 0;
    fgets(buffer, MAX_CHARS, infile);
    ninputs = sscanf(buffer, READ_FORMAT, range+i, ext+i,
		     radius+i, ext_bscat_ratio+i, ext_air+i,
		     ssa+i, g+i, ssa_air+i,
		     droplet_fraction+i, pristine_ice_fraction+i,
		     &nchars);
    if (ninputs < 4) {
      fprintf(stderr, "Error: garbled input at range gate %d: "
	      "only %d inputs read successfully\n", i+1, ninputs);
      exit(MS_INPUT_FILE_ERROR);
    }
    else if (ninputs < 10) {
      /* Characteristics of the phase function near backscatter have
	 not been entered: use defaults */
      droplet_fraction[i] = 0.0; /* Default: isotropic at 180deg */
      pristine_ice_fraction[i] = 0.0; /* Isotropic at 180deg */
      use_isotropic_pp = 1;
      if (ninputs < 8) {
	/* Wide-angle scattering properties have not been provided:
	   use small-angle multiple scattering only */
	config.wide_angle_algorithm = MS_WIDE_ANGLE_NONE;
	if (ninputs < 5) {
	  /* No air extinction: assume vacuum */
	  ext_air[i] = 0.0;
	  use_air_ext = 0;
	}
      }
    }
    else if (output_adjoint) {
      /* Read in the backscatter adjoint inputs */
      char* curbuf = buffer+nchars;
      char* newcurbuf = curbuf;
      int ifov;
      for (ifov = 0; ifov < instrument.nfov; ifov++) {
	bscat_AD[ifov*m+i] = strtod(curbuf, &newcurbuf);
	if (newcurbuf == curbuf) {
	  fprintf(stderr, "Error reading backscatter adjoint of field-of-view %d at range gate %d (should be value %d on the line)\n",
		  ifov+1, i+1, 11+ifov);
	  exit(MS_INPUT_FILE_ERROR);
	}
	else {
	  curbuf = newcurbuf;
	}
      }

      if (separate_bscat_air) {
	for (ifov = 0; ifov < instrument.nfov; ifov++) {
	  bscat_air_AD[ifov*m+i] = strtod(curbuf, &newcurbuf);
	  if (newcurbuf == curbuf) {
	    fprintf(stderr, "Error reading adjoint of air backscatter required with options \"-hsrl -adjoint\"\n");
	    exit(MS_INPUT_FILE_ERROR);
	  }
	  else {
	    curbuf = newcurbuf;
	  }
	}
      }
    }
  }

  if (automatically_configure) {
    if (instrument.wavelength > MS_RADAR_LIDAR_TRANSITION_WAVELENGTH) {
      /* Assume we have a radar */
      if (!manually_select_algorithms) {
	config.small_angle_algorithm = MS_SINGLE_SCATTERING;
      }
      if ((!manually_select_algorithms) || automatically_configure_tdts) {
	config.wide_angle_algorithm = MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE;
      }
      if (!manually_select_cbh) {
	config.coherent_backscatter_enhancement = 2.0;
      }
      instrument.receiver_type = MS_GAUSSIAN;
    }
    else {
      /* Assume we have a lidar */
    }
  }
  else if (automatically_configure_tdts) {
    if (instrument.wavelength > MS_RADAR_LIDAR_TRANSITION_WAVELENGTH) {
      config.wide_angle_algorithm = MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE;
    }
  }

  /* Allocate space for output */
  bscat = (ms_real*)malloc(sizeof(ms_real)*instrument.nfov*m);
  if (separate_bscat_air) {
    bscat_air = (ms_real*)malloc(sizeof(ms_real)*instrument.nfov*m);
  }
  if (!bscat || ((!bscat_air) && separate_bscat_air) ) {
    fprintf(stderr, "Error allocating space for apparent backscatter output\n");
    exit(MS_MEMORY_ALLOCATION_ERROR);
  }


  /* REPORT WHAT CALCULATIONS ARE BEING PERFORMED */
  if (!(config.options & MS_QUIET)) {
    ms_print_algorithms(config, instrument, range, use_isotropic_pp,
			separate_bscat_air, output_adjoint, output_jacobian,
			calc_jacobian, stderr);
  }

  if (m < n) {
    //    fprintf(stderr, "Warning: currently the number of samples (%d) cannot be fewer than\n"
    //	    "   the number of data points; setting them to be equal\n");*/
    n = m;
  }

  if (instrument.receiver_type == MS_GAUSSIAN
      && (config.small_angle_algorithm == MS_SMALL_ANGLE_PVC_ORIGINAL
	  || config.small_angle_algorithm == MS_SMALL_ANGLE_PVC_EXPLICIT)) {
    fprintf(stderr, "Warning: PVC algorithm will use a top-hat receiver pattern\n");
  }
  if (instrument.wavelength > MS_RADAR_LIDAR_TRANSITION_WAVELENGTH 
      && config.small_angle_algorithm > MS_SINGLE_SCATTERING) {
    fprintf(stderr, "Warning: wavelength greater than %g microns yet using small-angle scattering;\n"
	    "  should the \"-auto\" option be specified?\n",
	    MS_RADAR_LIDAR_TRANSITION_WAVELENGTH*1.0e6);
  }
  if (instrument.wavelength > MS_RADAR_LIDAR_TRANSITION_WAVELENGTH 
      && config.wide_angle_algorithm == MS_WIDE_ANGLE_TDTS_FORWARD_LOBE) {
    fprintf(stderr, "Warning: wavelength greater than %g microns yet using wide-angle scattering\n"
	    "  assuming a forward lobe; should the \"-auto\" option be specified?\n",
	    	    MS_RADAR_LIDAR_TRANSITION_WAVELENGTH*1.0e6);
  }


  /* CHECK RANGE-GATE SPACING */
  if ((!ms_range_spacing_is_regular(n, range, 
				   MS_RANGE_SPACING_TOLERANCE))
      && config.wide_angle_algorithm != MS_WIDE_ANGLE_NONE) {
    fprintf(stderr, 
	    "  The range-gate spacing is not constant to within a tolerance of 5%%,\n"
	    "  so will be interpolated on to a regular grid for the purposes of wide-angle scattering\n");
  }

  /* RUN ALGORITHM */
//  if (config.options & MS_CRUDE_INTEGRATION) {
//    ext_air = ssa_air = NULL;
//  }
  if (!output_jacobian) {
    /* Call multiscatter algorithm */
    for (i = 0; i < nrepeats; i++) {
      if (!output_adjoint) {
	CHECK(multiscatter(n, m, &config, instrument, surface, range, radius,
			   ext, ssa, g, ext_bscat_ratio, ext_air, ssa_air,
			   droplet_fraction, pristine_ice_fraction,
			   bscat, bscat_air));
      }
      else {
	//	ms_real radius_AD[n];
	ms_real ext_air_AD[n];
	ms_real droplet_fraction_AD[n];
	ms_real pristine_ice_fraction_AD[n];
	int i;
	for (i = 0; i < n; i++) {
	  ext_air_AD[i] = droplet_fraction_AD[i] = pristine_ice_fraction_AD[i]
	    = 0.0;
	}
	CHECK(multiscatter_AD(n, m, &config, instrument, surface, range, radius,
			      ext, ssa, g, ext_bscat_ratio, ext_air, ssa_air,
			      droplet_fraction, pristine_ice_fraction,
			      bscat, bscat_air,
			      bscat_AD, bscat_air_AD,
			      radius_AD, ext_AD, ssa_AD, g_AD, 
			      ext_bscat_ratio_AD, ext_air_AD, 
			      droplet_fraction_AD, pristine_ice_fraction_AD));
      }
      config.options |= MS_QUIET;	  
    }

    /* Output results */
    if (config.options & MS_OUTPUT_DISTRIBUTION) {
      ms_print_distribution(config, stdout);
    }
    else {
      for (i = 0; i < n; i++) {
	int ifov;
	/* Print the basic variables */
	fprintf(stdout, "%d %g %g %g %14.9g",
		i+1, 
		range[i],
		ext[i],
		radius[i],
		bscat[i]);
	if (instrument.nfov > 1) {
	  /* More than one field-of-view: loop over the others and
	     print the backscatters */
	  for (ifov = 1; ifov < instrument.nfov; ifov++) {
	    fprintf(stdout, " %14.9g", bscat[i + ifov*m]);
	  }
	}

	/* Print the air backscatter */
	if (separate_bscat_air) {
	  for (ifov = 0; ifov < instrument.nfov; ifov++) {
	    fprintf(stdout, " %14.9g", bscat_air[i + ifov*m]);
	  }
	}
	
	/* Print the mean path lag (experimental) */
	if (config.small_angle_algorithm == MS_SMALL_ANGLE_PVC_FAST_LAG) {
	  fprintf(stdout, " %g", config.small_angle_lag[i]);
	}

	/* Print the adjoint */
	if (output_adjoint) {
	  fprintf(stdout, " %14.9g %14.9g %14.9g %14.9g %14.9g",
		  ext_AD[i], ssa_AD[i], g_AD[i], ext_bscat_ratio_AD[i], radius_AD[i]);
	}
	/*
	if ((config.options & MS_OUTPUT_BSCATS) && (config.options & MS_EXPLICIT_QSA)) {
	  fprintf(stdout, " %g %g %g 0.0 0.0 0.0",
		  ms_bscat_single[i],
		  ms_bscat_double[i],
		  ms_bscat_multi[i]);
	}
	*/
	fprintf(stdout, "\n");
      }
      if (m > n) {
	ms_real drange = range[1]-range[0];
	for (i = n; i < m; i++) {
	  fprintf(stdout, "%d %g %g %g %18.9g",
		  i+1, 
		  range[n-1]+drange*(i-n+1),
		  0.0,
		  0.0,
		  bscat[i]);
	  if (instrument.nfov > 1) {
	    /* More than one field-of-view: loop over the others and
	       print the backscatters */
	    int ifov;
	    for (ifov = 1; ifov < instrument.nfov; ifov++) {
	      fprintf(stdout, " %14.9g", bscat[i + ifov*m]);
	    }
	  }
	  
	  /* Print the air backscatter (zero) */
	  if (separate_bscat_air) {
	    fprintf(stdout, " 0.0");
	  }
	  if (output_adjoint) {
	    fprintf(stdout, " 0.0 0.0 0.0 0.0");
	  }
	  fprintf(stdout, "\n");
	}
      }
    }
  }

  /* JACOBIAN CALCULATION */

  if (calc_jacobian) {
    ms_real *jacobian_data
      = (ms_real*)malloc(4*m*instrument.nfov*n*sizeof(ms_real));
    ms_real *jacobian_air_data = NULL;
    ms_real **d_bscat_d_ext = (ms_real**)malloc(n*sizeof(ms_real*));
    ms_real **d_bscat_d_ssa = (ms_real**)malloc(n*sizeof(ms_real*));
    ms_real **d_bscat_d_g = (ms_real**)malloc(n*sizeof(ms_real*));
    ms_real **d_bscat_d_radius = (ms_real**)malloc(n*sizeof(ms_real*));
    ms_real *d_bscat_d_ext_bscat_ratio
      = (ms_real*)malloc(m*instrument.nfov*m*sizeof(ms_real));
    ms_real **d_bscat_air_d_ext = NULL;
    ms_real **d_bscat_air_d_ssa = NULL;
    ms_real **d_bscat_air_d_g = NULL;
    ms_real **d_bscat_air_d_radius = NULL;
    int *i_xy = (int*)malloc(m*instrument.nfov*sizeof(int));
    for (i = 0; i < m*instrument.nfov; i++) {
      i_xy[i] = i;
    }
    for (i = 0; i < n; i++) {
      d_bscat_d_ext[i] =  jacobian_data + i*m*instrument.nfov;
      d_bscat_d_ssa[i] =  jacobian_data + (i+n)*m*instrument.nfov;
      d_bscat_d_g[i]   =  jacobian_data + (i+2*n)*m*instrument.nfov;
      d_bscat_d_radius[i]=jacobian_data + (i+3*n)*m*instrument.nfov;
    }

    if (bscat_air) {
      jacobian_air_data 
	= (ms_real*)malloc(4*m*instrument.nfov*n*sizeof(ms_real));
      d_bscat_air_d_ext = (ms_real**)malloc(n*sizeof(ms_real*));
      d_bscat_air_d_ssa = (ms_real**)malloc(n*sizeof(ms_real*));
      d_bscat_air_d_g = (ms_real**)malloc(n*sizeof(ms_real*));
      d_bscat_air_d_radius = (ms_real**)malloc(n*sizeof(ms_real*));
      for (i = 0; i < n; i++) {
	d_bscat_air_d_ext[i] =  jacobian_air_data + i*m*instrument.nfov;
	d_bscat_air_d_ssa[i] =  jacobian_air_data + (i+n)*m*instrument.nfov;
	d_bscat_air_d_g[i]   =  jacobian_air_data + (i+2*n)*m*instrument.nfov;
	d_bscat_air_d_radius[i]=jacobian_air_data + (i+3*n)*m*instrument.nfov;
      }
    }
    for (i = 0; i < nrepeats; i++) {
      if (calc_jacobian == 1) {
	CHECK(ms_jacobian_linear(n, m, 
		 &config, instrument, surface,
		 range, radius, ext, ssa, g,
		 ext_bscat_ratio, ext_air, ssa_air,
		 droplet_fraction, pristine_ice_fraction,
		 n, m*instrument.nfov, i_xy, i_xy, 
		 bscat, bscat_air,
		 d_bscat_d_ext, d_bscat_d_ssa, d_bscat_d_g,
		 d_bscat_d_radius, d_bscat_d_ext_bscat_ratio,
		 d_bscat_air_d_ext, d_bscat_air_d_ssa, 
		 d_bscat_air_d_g, d_bscat_air_d_radius));
      }
      else {
	/* Only require Jacobian with respect to extinction */
	CHECK(ms_jacobian_linear(n, m, 
		   &config, instrument, surface,
		   range, radius, ext, ssa, g,
		   ext_bscat_ratio, ext_air, ssa_air,
		   droplet_fraction, pristine_ice_fraction,
		   n, m*instrument.nfov, i_xy, i_xy, 
		   bscat, bscat_air,
		   d_bscat_d_ext, NULL, NULL, NULL, NULL,
		   d_bscat_air_d_ext, NULL, NULL, NULL));
      }
    }	    

    if (output_adjoint == 2) {
      /* Check the adjoint against the numerical Jacobian */
      ms_real ext_AD_num[n];
      ms_real ssa_AD_num[n];
      ms_real radius_AD_num[n];
      ms_real g_AD_num[n];
      ms_real ext_bscat_ratio_AD_num[n];
      int i;
      for (i = 0; i < n; i++) {
	ext_AD_num[i] = ssa_AD_num[i]
	  = g_AD_num[i] = ext_bscat_ratio_AD_num[i] 
	  = radius_AD_num[i] = 0.0;
      }

      for (i = 0; i < m*instrument.nfov; i++) {
	if (bscat_AD[i] != 0.0) {
	  int j;
	  for (j = 0; j < n; j++) {
	    ext_AD_num[j] += d_bscat_d_ext[j][i]*bscat_AD[i];
	  }
	}
	if (calc_jacobian == 1) {
	  /* Jacobian is available with respect to other variables */
	  int j;
	  for (j = 0; j < n; j++) {
	      ssa_AD_num[j] += d_bscat_d_ssa[j][i]*bscat_AD[i];
	      g_AD_num[j]   += d_bscat_d_g[j][i]*bscat_AD[i];
	      radius_AD_num[j] += d_bscat_d_radius[j][i]*bscat_AD[i];
	  }
	  j = i % m;
	  if (j < n) {
	    ext_bscat_ratio_AD_num[j]
	      += d_bscat_d_ext_bscat_ratio[j]*bscat_AD[i];
	  }
	}
      }
      
      if (bscat_air) {
	for (i = 0; i < m*instrument.nfov; i++) {
	  if (bscat_air_AD[i] != 0.0) {
	    int j;
	    for (j = 0; j < n; j++) {
	      ext_AD_num[j] += d_bscat_air_d_ext[j][i]*bscat_air_AD[i];
	    }
	    if (calc_jacobian == 1) {
	      /* Jacobian is available with respect to other variables */
	      for (j = 0; j < n; j++) {
		ssa_AD_num[j] += d_bscat_air_d_ssa[j][i]*bscat_air_AD[i];
		g_AD_num[j]   += d_bscat_air_d_g[j][i]*bscat_air_AD[i];
		radius_AD_num[j] += d_bscat_air_d_radius[j][i]*bscat_air_AD[i];
	      }
	    }
	  }
	}
      }

      fprintf(stderr, "Result of adjoint evaluation using a numerical Jacobian:\n");
      if (calc_jacobian == 1) {
	fprintf(stderr, " ext_AD    ssa_AD      g_AD ext_bscat_ratio_AD radius_AD bscat_AD(%d FOVs) bscat_air_AD(%d FOVs)\n",
	      instrument.nfov, instrument.nfov);
      }
      else {
	fprintf(stderr, " ext_AD  bscat_AD(%d FOVs) bscat_air_AD(%d FOVs)\n",
	      instrument.nfov, instrument.nfov);
      }
      for (i = 0; i < n; i++) {
	int ifov;
	print_error_AD(stdout, ext_AD_num[i], ext_AD[i]);
	if (calc_jacobian == 1) {
	  print_error_AD(stdout, ssa_AD_num[i], ssa_AD[i]);
	  print_error_AD(stdout, g_AD_num[i], g_AD[i]);
	  print_error_AD(stdout, ext_bscat_ratio_AD_num[i], 
			 ext_bscat_ratio_AD[i]);
	  print_error_AD(stdout, radius_AD_num[i], radius_AD[i]);
	}
	fprintf(stdout, "        ");
	for (ifov = 0; ifov < instrument.nfov; ifov++) {
	  fprintf(stdout, "%g ", bscat_AD[i+ifov*m]);
	}
	fprintf(stdout, "        ");
	if (bscat_air_AD) {
	  for (ifov = 0; ifov < instrument.nfov; ifov++) {
	    fprintf(stdout, "%g ", bscat_air_AD[i+ifov*m]);
	  }
	}
	else {
	  fprintf(stdout, "N/A  ");
	}
	fprintf(stdout, "\n");
      }
    }

    if (output_jacobian) {
      /* First Jacobian with respect to extinction */
      int j;
      for (j = 0; j < n; j++) {
	int i;
	for (i = 0; i < m*instrument.nfov; i++) {
	  fprintf(stdout, " %g", d_bscat_d_ext[j][i]);
	}
	if (bscat_air) {
	  for (i = 0; i < m*instrument.nfov; i++) {
	    fprintf(stdout, " %g", d_bscat_air_d_ext[j][i]);
	  }
	}
	fprintf(stdout, "\n");
      }
      if (calc_jacobian == 1) {
	/* Jacobian is available with respect to other variables */
	/* Jacobian with respect to single scattering albedo */
	for (j = 0; j < n; j++) {
	  int i;
	  for (i = 0; i < m*instrument.nfov; i++) {
	    fprintf(stdout, " %g", d_bscat_d_ssa[j][i]);
	  }
	  if (bscat_air) {
	    for (i = 0; i < m*instrument.nfov; i++) {
	      fprintf(stdout, " %g", d_bscat_air_d_ssa[j][i]);
	    }
	  }
	  fprintf(stdout, "\n");
	}
	/* Jacobian with respect to asymmetry factor */
	for (j = 0; j < n; j++) {
	  int i;
	  for (i = 0; i < m*instrument.nfov; i++) {
	    fprintf(stdout, " %g", d_bscat_d_g[j][i]);
	  }
	  if (bscat_air) {
	    for (i = 0; i < m*instrument.nfov; i++) {
	      fprintf(stdout, " %g", d_bscat_air_d_g[j][i]);
	    }
	  }
	  fprintf(stdout, "\n");
	}
	/* Jacobian with respect to radius */
	for (j = 0; j < n; j++) {
	  int i;
	  for (i = 0; i < m*instrument.nfov; i++) {
	    fprintf(stdout, " %g", d_bscat_d_radius[j][i]);
	  }
	  if (bscat_air) {
	    for (i = 0; i < m*instrument.nfov; i++) {
	      fprintf(stdout, " %g", d_bscat_air_d_radius[j][i]);
	    }
	  }
	  fprintf(stdout, "\n");
	}
	/* Jacobian with respect to extinction-to-backscatter ratio */
	for (i = 0; i < m*instrument.nfov; i++) {
	  if (i < n) {
	    fprintf(stdout, " %g", d_bscat_d_ext_bscat_ratio[i]);
	  }
	  else {
	    fprintf(stdout, " 0");
	  }
	  if (bscat_air) {
	    for (i = 0; i < m*instrument.nfov; i++) {
	      fprintf(stdout, " 0");
	    }
	  }
	}
	fprintf(stdout, "\n");
      }
    }
  }

  if (0) {
    /* Calculate the Jacobian */
    ms_real *jacobian_data = (ms_real*)malloc((2*n+1)*n*sizeof(ms_real));
    ms_real **d_ln_beta_d_ln_ext = (ms_real**)malloc(n*sizeof(ms_real *));
    ms_real **d_ln_beta_d_ln_radius = (ms_real**)malloc(n*sizeof(ms_real *));
    ms_real *d_ln_beta_d_ln_ext_bscat_ratio = jacobian_data + 2*n*n;
    int *i_xy = (int*)malloc(n*sizeof(int));
    int i, j;
    for (i = 0; i < n; i++) {
      d_ln_beta_d_ln_ext[i] = jacobian_data + i*n;
      d_ln_beta_d_ln_radius[i] = jacobian_data + (n+i)*n;
      i_xy[i] = i;
    }
    for (i = 0; i < nrepeats; i++) {
      CHECK(ms_jacobian(n, m, &config, instrument, surface,
				  range, radius, ext, ssa, g,
				  ext_bscat_ratio, ext_air, ssa_air,
				  droplet_fraction, pristine_ice_fraction,
				  n, n, i_xy, i_xy, 
				  bscat,
				  d_ln_beta_d_ln_ext, d_ln_beta_d_ln_radius, 
				  d_ln_beta_d_ln_ext_bscat_ratio));
    }
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
	fprintf(stdout, " %g", d_ln_beta_d_ln_ext[i][j]);
      }
      fprintf(stdout, "\n");
    }
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
	fprintf(stdout, " %g", d_ln_beta_d_ln_radius[i][j]);
      }
      fprintf(stdout, "\n");
    }
    for (i = 0; i < n; i++) {
      fprintf(stdout, " %g", d_ln_beta_d_ln_ext_bscat_ratio[i]);
    }
    fprintf(stdout, "\n");
  }

  if (print_stats) {
    ms_print_stats(config, stderr);
  }

  exit(MS_SUCCESS);
}
