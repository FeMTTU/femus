#ifndef __mgqrule_h__
#define __mgqrule_h__


#include <string>
#include <vector>

#include "Typedefs.hpp"
#include "VBTypeEnum.hpp"
#include "GeomEl.hpp"

namespace femus {



class QRule  {

  public:

//GeomEl =========
  GeomEl _geomel;

//Quadrature =========
    std::string  _qrule_type;
    uint         _NoGaussVB;
    std::vector<double>   _weightVB;
    
    
     QRule(GeomEl geomel_in);
    
  
    
};


} //end namespace femus



#endif