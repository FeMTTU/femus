#ifndef _typedefs_cfg_h_
#define _typedefs_cfg_h_


//this is mostly for primitives data types


typedef unsigned int uint;   //needed everywhere

#include "FemusExtLib_conf.hpp"

//******************** 
 #ifdef LM_REAL
 #include "libmesh/libmesh_common.h"
 #else
 typedef double Real;
 #endif
//******************



  
#endif