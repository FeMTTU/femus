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


enum FEFamiliesAsAList {
    CONTINUOUS_LINEAR      = 0,
    CONTINUOUS_SERENDIPITY,
    CONTINUOUS_BIQUADRATIC,
    DISCONTINUOUS_CONSTANT,
    DISCONTINUOUS_LINEAR
};



#endif
