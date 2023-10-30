#ifndef __femus_enums_FETypeEnum_hpp__
#define __femus_enums_FETypeEnum_hpp__

#include <string>
#include <vector>

namespace femus {

  static const std::vector< std::string >  fe_fams    = {"linear", "quadratic", "biquadratic", "constant", "disc_linear"};

}

#define NFE_FAMS 5
#define NFE_FAMS_C0_LAGRANGE  3  

#define BIQUADR_FE   2
#define LINEAR_FE    0



enum  FEType {QQ=0, LL, KK };   ///@deprecated

#define QL_NODES 2  ///@deprecated
#define QL 3        ///@deprecated

#define  MESH_MAPPING_FE  1  ///@deprecated linear 

#define MESH_ORDER 0        ///@deprecated  biq


#endif
