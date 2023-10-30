#ifndef __femus_enums_FETypeEnum_list_hpp__
#define __femus_enums_FETypeEnum_list_hpp__


// This file treats the FE families in the library as a one-indexed list


#include <string>
#include <vector>

namespace femus {

  static const std::vector< std::string >  fe_fams    = {"linear", "quadratic", "biquadratic", "constant", "disc_linear"};

}

#define NFE_FAMS 5
#define NFE_FAMS_C_ZERO_LAGRANGE  3  

#define LINEAR_FE    0
#define BIQUADR_FE   2





#endif
