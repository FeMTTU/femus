#ifndef __femus_enums_FETypeEnum_hpp__
#define __femus_enums_FETypeEnum_hpp__

#include <string>
#include <vector>

namespace femus {

  static const std::vector< std::string >  fe_fams    = {"linear", "quadratic", "biquadratic", "constant", "disc_linear"};

}

#define NFE_FAMS 5



enum  FEType {QQ=0, LL, KK };

#define QL_NODES 2
#define QL 3

#define BIQUADR_FE   2
#define LINEAR_FE    0


#define  MESH_MAPPING_FE  1  //linear
#define MESH_ORDER 0          //biq

#endif
