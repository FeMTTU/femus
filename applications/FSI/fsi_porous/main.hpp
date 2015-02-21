#ifndef __main_hpp__
#define __main_hpp__


#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<cstring>

#include <vector>
#include <map>

#include "petscksp.h"
#include "petscvec.h" 


using std::vector;
using std::map;
using std::cout;
using std::endl;

///////////////// application specific includes

//compile time
#define SPACEDIM 3

#define MAX_EL_NODES   27 //3D, biquadratic, not more than that


#define N_QTIES  3
#define IDX_P   0
#define IDX_D   1 
#define IDX_VEL 2


#define NVAR_D 3
#define NVAR_VEL 3
#define NVAR_P 1


// NLEVS is runtime now!



#endif

//YOUNG_FRAC
//  1e4: it doesnt cross
//
// 0.8e4: it crosses
// 0.4e4: it crosses SORT OF TANGENTIAL
// 0.2e4: it crosses

// E_frac = 0.5e4 vs 1.e10: Productivity index significantly different