#ifndef __feelembase_h__
#define __feelembase_h__

#include <string>
#include <vector>

#include "Typedefs.hpp"
#include "VBTypeEnum.hpp"
#include "GeomEl.hpp"
#include "QRule.hpp"
#include "ElemType.hpp"

namespace femus {



//Basic finite element for VB
//it is not quadratic and linear altogether, it is only quadratic or linear
// The CONNECTION between GEOM ELEMENT, FE And QUADRATURE is to be cleared out
// The generation of the sparsity pattern in particular starts with mesh - geom-element - feelem (nodetodof)
// one understands that for instance if you have a Quad8 mesh and you want a Quad9 finite element

class FEElemBase  {

public:

    FEElemBase(std::vector<GeomEl> geomel_in);
    virtual ~FEElemBase();

// GeomEl ======
    std::vector<GeomEl>  _geomel;   //VB

// FE ==========
    static  FEElemBase* build(std::vector<GeomEl> geomel_in, const uint order);
    
// Quadrature ==
    std::vector<QRule> _qrule;   //VB
    void AssociateQRule(std::vector<QRule> qrule_in);
    double**      _phi_mapVBGD[VB];
    double** _dphidxez_mapVBGD[VB];
    void evaluate_shape_at_qp(const uint order);  /*This argument will be removed*/
    
// Multigrid ======
    virtual float get_embedding_matrix(const uint,const uint,const uint) = 0;
    virtual double get_prol(const uint) = 0;

// ====================
  std::vector<elem_type*> _myelems;    //VB
  typedef double* (elem_type::*_FunctionPointerTwo)(const unsigned & ig) const; //declaring the FunctionPointer type
  std::vector< std::vector<_FunctionPointerTwo> > _DphiptrTwo;
    
    

};


//the fact of defining CHILDREN of an element is more a GEOMETRIC THING
//in the sense that it can just depend on the geometry... but having 
//children of Hex27 or children of Hex20 is not the same for the FE dofs
//also the number of children can also be related to the REFINEMENT ALGORITHM.
//the refinement is associated to the MESH rather than to the DOFS i think
//we use midpoint refinement.
//this must be static, so that you can call it even without an instantiation!

//what is the difference between a STATIC function within a class 
//and a static function outside?
//in this way the class is like a "namespace" for the function...


} //end namespace femus



#endif