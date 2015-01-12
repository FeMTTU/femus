#ifndef __mggeomel_h__
#define __mggeomel_h__


#include <string>

#include "Typedefs.hpp"
#include "FETypeEnum.hpp"


namespace femus {

// Merge this in Elem... but wait, this is all ABSTRACT, what about Elem?

class GeomEl  {

public:

     GeomEl(const std::string geomel_id, const uint mesh_order);
    ~GeomEl();
    
    uint _dim;               /*THIS CANNOT BE CONST OTHERWISE I CANNOT DO A VECTOR WITH A PUSH_BACK (or I should create a copy constructor...)*/
    std::string _geomel_id;      
    uint      _elnds;        // number of nodes of the element
    std::string name;        // element name: quadratic and corresponding linear 

//===== Multigrid   
    uint n_se;                    ///< number of subelements 
    
    //from linear to quadratic
//     static const double _Prol[NNDS*NNDSL];  //actually NNDSL should be NDOF_P: one part is geometric, one part is mathematic


};



} //end namespace femus


#endif