#ifndef __currelem_h__
#define __currelem_h__


#include "VBType_enum.hpp"
#include "FEType_enum.hpp"

#include "DenseVector.hpp"
#include "DenseMatrix.hpp"

class EqnBase;
class EquationsMap;
class QuantityLocal;


// Ok, in this case if we want to reach all the other classes we just 
// use the eq pointer, through which we have the eqn map, through which 
// we have everything...
//Through an equation we get to the equationsmap which gives us access to everything we need
//Ok, I need to pass the Equation, and also I will pass the ref to EqMap, 
// which is protected in the EqnBase class and I want to leave it protected

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
 //TODO here we only need the GEOMETRIC PART of the CurrElem, so maybe we will split 
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

//So, so far we have the CurrElem class is 
//split into a CURR GEOMETRIC and a CURR EQUATION part.
//The curr geometric is basically filled with the MESH class
//The curr fe is basically filled with the EQUATION class

//can we make a currelem that is not based on the mesh?
//Let us now just make this one

  class CurrElem {

  public:
    
    CurrElem(EqnBase&,EquationsMap& e_map_in);
   ~CurrElem();

   EqnBase & _eqn;  //con questo puoi accedere a dati e funzioni DEL PADRE, NON al FIGLIO
   EquationsMap & _eqnmap;

// ========================================================================================
//========== ELEMENT: Current Geometric Element (SERVICE)  ========================
     uint  **_el_conn;             /// vector of the global nodes for that element         [VB][NNDS];
     uint    _vol_iel_DofObj[VB];  /// i need to put the element also. both VV and BB
   double  **_xx_nds;              /// vector of the node coordinates for that element     [VB][_spacedimension*NNDS];
   double  **_el_xm;               /// element center point                               [VB][_spacedimension];

    //==== REAL ELEMENT properties: connectivity, coordinates, ...  need the MESH BASICALLY
    void  get_el_nod_conn_lev_subd(const uint vb,const uint Level,const uint isubd_in,const uint iel) const;
    void  get_el_DofObj_lev_subd(const uint vb,const uint Level,const uint isubd_in,const uint iel);  //TODO this is not const because _vol_iel_DofObj is NOT A POINTER! 
    void  get_el_ctr(const uint bdry) const;                                                          //TODO notice that this is not changing the POINTER , so it is const!
    void  get_el_orient(const uint vb) const;
    void  ConvertElemCoordsToMappingOrd(const uint vb,QuantityLocal& myvect) const;
   
// ========================================================================================
//========== ELEMENT: Current "EQUATION" Element (ql are TOGETHER ) ========================               
  uint                   _el_n_dofs[VB];
  std::vector<uint> _el_dof_indices[VB];
  uint*                  _bc_eldofs[VB]; //So the element must be aware of the BC of the equation
  DenseMatrix                 _KeM[VB]; 
  DenseVector                 _FeM[VB];
  
 void GetElDofsBc(const uint vbfl, const uint Level);  //needs the EQUATION basically
 
   
  };
  

#endif