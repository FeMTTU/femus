#ifndef __mggeomel_h__
#define __mggeomel_h__


#include <string>

#include "Typedefs.hpp"
#include "VBTypeEnum.hpp"
#include "FETypeEnum.hpp"


namespace femus {



class GeomEl  {

public:

     GeomEl(const uint dim_in,const uint geomel_type);
    ~GeomEl();
    
    uint _geomel_type;   /*THIS CANNOT BE CONST OTHERWISE I CANNOT DO A VECTOR WITH A PUSH_BACK (or I should create a copy constructor...)*/
    uint _dim;           /*THIS CANNOT BE CONST OTHERWISE I CANNOT DO A VECTOR WITH A PUSH_BACK (or I should create a copy constructor...)*/
    uint _elnds[QL_NODES]; // I SWITCHED THIS ONE: VB first, QL next!!!   //number of nodes of one element [QL]
    std::string name;             ///< element name 
    std::string pname;            ///< print element name for XDMF print (linear)

//===== Multigrid   
    uint n_se;                    ///< number of subelements 
    
    //from linear to quadratic
//     static const double _Prol[NNDS*NNDSL];  //actually NNDSL should be NDOF_P: one part is geometric, one part is mathematic


};



} //end namespace femus


#endif