#ifndef __femus_meshGencase_FEElemBase_hpp__
#define __femus_meshGencase_FEElemBase_hpp__

#include <string>
#include <vector>

#include "Typedefs.hpp"
#include "VBTypeEnum.hpp"
#include "GaussPoints.hpp"
#include "ElemType.hpp"

namespace femus {



class FEElemBase  {

public:

    FEElemBase();
    virtual ~FEElemBase();

    static  FEElemBase* build(const std::string geomel_id_in, const uint fe_family);

    virtual unsigned int get_dimension() const { std::cout << "Not implemented FE" << std::endl; abort(); };
    virtual unsigned int n_nodes()       const { std::cout << "Not implemented FE" << std::endl; abort(); };
    virtual std::string  get_name_med()  const { std::cout << "Not implemented FE" << std::endl; abort(); };
    virtual std::string  get_name_xdmf() const { std::cout << "Not implemented FE" << std::endl; abort(); };
    
// Multigrid ===   //only VV
    virtual float get_embedding_matrix(const uint,const uint,const uint) = 0;  //should be geom_el
    virtual double get_prol(const uint) = 0;
 
    
};


} //end namespace femus



#endif
