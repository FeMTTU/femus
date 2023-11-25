#ifndef __femus_enums_FElemTypeEnum_hpp__
#define __femus_enums_FElemTypeEnum_hpp__


#include <string>
#include <vector>


enum  FEFamily {
    LAGRANGE = 0,
    DISCONTINUOUS_POLYNOMIAL,
    WEAK_GALERKIN
};

constexpr unsigned int FEFamily_count = WEAK_GALERKIN - LAGRANGE + 1;


enum FEOrder {
    ZERO = 0,
    FIRST,
    SERENDIPITY,
    SECOND
};



namespace femus {

  static const std::vector< std::string >  fe_families          = {"Continuous Lagrange", "Discontinuous Polynomial"};
  
}



#endif
