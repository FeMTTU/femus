/*=========================================================================

 Program: FEMUS
 Module: CurrentElem
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __currelem_h__
#define __currelem_h__


#include "FETypeEnum.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "MeshTwo.hpp"
#include "DenseVector.hpp"
#include "DenseMatrix.hpp"


namespace femus {



class SystemTwo;
class MultiLevelProblemTwo;
class CurrentQuantity;



  class CurrentElem {

  public:
    
    CurrentElem(const uint vb, const SystemTwo*, const MeshTwo& mesh, const std::vector< std::vector<elem_type*> >  & elem_type);
   ~CurrentElem();

    inline const uint  GetVb() const {
      return _mesh.get_dim() - _dim;
    }
    
    inline const uint  GetDim() const {
      return _dim;
    }
    
    inline const double*  GetMidpoint() const {
      return _el_xm;
    }
    
    inline const uint*  GetConn() const {
      return _el_conn;
    }
    
    inline uint  GetVolIel() const {
      return _vol_iel_DofObj;
    }
    
    inline const double*  GetNodeCoords() const {
      return _xx_nds;
    }
    
    inline std::vector<uint>  GetDofIndices() const {
      return _el_dof_indices;
    }
    
    //TODO here i have to return non-const because of a function that i will change...
    inline uint*  GetBCDofFlag() const {
      return _bc_eldofs;
    }
    
    /** Returns a reference to the element rhs */
    inline DenseVector & Rhs() {
      return _FeM;
    }
    
    /** Returns a reference to the element matrix */
    inline DenseMatrix & Mat() {
      return _KeM;
    }
    
    void  set_el_nod_conn_lev_subd(const uint Level,const uint isubd_in,const uint iel);
    
    //TODO notice that this is not changing the POINTER , so it is const!
    void  SetMidpoint() const;
    
    void  PrintOrientation() const;
    
    void  ConvertElemCoordsToMappingOrd(CurrentQuantity& myvect) const;
    
    /** needs the EQUATION basically */
    void  SetElDofsBc(const uint Level);
 
    //TODO make these private
//========== Equation-related ========================               
  const SystemTwo * _eqn;  //con questo puoi accedere a dati e funzioni DEL PADRE, NON al FIGLIO
  const MeshTwo & _mesh;
  const std::vector< std::vector<elem_type*> >  &  _elem_type;
  
  private:
    
// ========================================================================================
//========== Current "EQUATION" Element (ql are TOGETHER ): needs the EQUATION ========================               
  DenseMatrix                  _KeM; 
  DenseVector                  _FeM;
  uint                   _el_n_dofs;
  std::vector<uint> _el_dof_indices;
  uint*                  _bc_eldofs; //So the element must be aware of the BC of the equation
  
// ========================================================================================
//==========  Current Geometric Element:  needs the MESH  ========================
     uint  *_el_conn;             /// vector of the global nodes for that element         [NNDS];
     uint    _vol_iel_DofObj;     /// i need to put the element also.
   double  *_xx_nds;              /// vector of the node coordinates for that element     [_spacedimension*NNDS];
   double  *_el_xm;               /// element center point                                [_spacedimension];
   const uint _dim;         //spatial dimension of the current element (can be different from the mesh dimension!)
   const uint _mesh_vb;     //index for the mesh
    
  };
  


} //end namespace femus

// Ok, in this case if we want to reach all the other classes we just 
// use the eq pointer, through which we have the eqn map, through which 
// we have everything...
//Through an equation we get to the equationsmap which gives us access to everything we need
//Ok, I need to pass the Equation, and also I will pass the ref to EqMap, 
// which is protected in the SystemTwo class and I want to leave it protected

//The geometry does not depend on the specific equation
//The rest yes

// CurrentElement and CurrentGaussPoint are supposed to be "External containers"
// to hold the informations needed for those objects

//Since we structured the GenMatRhs in terms of VB, then the current elem and gauss point
//are going to have all the VB structure... maybe we could do an object for only VV
// and one for only BB, maybe templatized. So far we go ahead like this.

//here there should be two pointers:
// to the ABSTRACT GEometrical Element
// and to the ABSTRACT FE Element

//Basically earlier we were using the Mesh and the Equation to handle what was needed by the current element.
//Now we are using the current element and current gauss to call from the big classess,
// Mesh and Equation, what is needed for the assemblying

 //alright, here I need a Current Element, to perform the loop
 //the current element only needs the eqn map, so we can use it everywhere
 //TODO here we only need the GEOMETRIC PART of the CurrentElem, so maybe we will split 
 //between the CurrGeomEl and the CurrFEEl
//in fact this is used for many element loop, but just to retrieve the geometrical properties
//like coords, middle point, etc.

//Ok, so far we have VB, so if we want only BB or VV we should actually make a vector 
// of two functions like this
//TODO i should TEMPLATIZE over VB!!! YES IT IS TRUE, IN FACT BASICALLY EVERYTHING 
//DEPENDS ON VB, and when you use one you dont use the other!!!

//=======================
//ok here we would need a "REFERENCE REAL ELEMENT"
//the current element contains the absolute coordinates, 
//but not in the reference frame
//ok now i want to set the element center but NOT BASED ON THE MESH

//So, so far we have the CurrentElem class is 
//split into a CURR GEOMETRIC and a CURR EQUATION part.
//The curr geometric is basically filled with the MESH class
//The curr fe is basically filled with the EQUATION class

//can we make a currelem that is not based on the mesh?
//Let us now just make this one


#endif