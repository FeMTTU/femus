#ifndef __quantitylocal_h__
#define __quantitylocal_h__

#include <vector>
#include "VBTypeEnum.hpp"
#include "CurrGaussPoint.hpp"

// The main difference to consider is  Vect WITH    Quantity
//                                 and Vect WITHOUT Quantity 
// Once you have the Quantity, everything else comes from that association 
////////////////
 //we've seen that these things stay close to each other

 //in every equation there may be Vect with associated quantity or not.
 //So not necessarily you need to define a quantity to associate,
 //but it may be very helpful
 
 //*********
 // every Vect MAY OR MAY NOT HAVE A QUANTITY
 // A QUANTITY MAY OR MAY NOT HAVE AN EQUATION
 // IF A QUANTITY DOES NOT HAVE AN EQUATION, THEN IT MUST HAVE A FUNCTION
 
 //Since the Vect may or may not have a Quantity, 
 //I don't pass the quantity to the constructor, but later
 
 //Ok, I need the list of Abstract FE Elements.
//Now, Vect is first of all associated with a Quantity
//but wait, what if it is not associated to any Quantity?
//TODO then you need other ways to get to the FEMap.
//The quantity must have the FEMap...
//Then Vect must be used somewhere where the _FEMap can be reached
 
//TODO the Vect class should be called the VectGaussPoint, because it mostly holds GAUSS POINT VALUES ! 
 
 //Ok, so Vect can have a quantity with equation,
 //or a quantity without equation,
 //or no quantity at all
 
 //Vect has a LOCAL nature, in the sense that it has the ELEMENT dofs, and also the CURRENT GAUSS values
 //for several operators
 
 //TODO TODO TODO this class must be better redefined 
 
 
 class EqnBase;
 class Quantity;
 class CurrGaussPointBase;
 class CurrElem;

  class QuantityLocal {
    
  public:
    
     QuantityLocal(CurrGaussPointBase &, CurrElem &);
    ~QuantityLocal();

    
    double*  _val_g;                                      //NEED TO ALLOCATE THIS ONE BEFORE IF YOU USE IT  //WHY ARE WE NOT ALLOCATING IN THE CONSTRUCTOR?!?
    double*  _val_g3D;   //for cross products             //NEED TO ALLOCATE THIS ONE BEFORE IF YOU USE IT
    double*  _val_dofs;   //NEED TO ALLOCATE THIS ONE BEFORE IF YOU USE IT
    double*  _val_dofs3D;   //NEED TO ALLOCATE THIS ONE BEFORE IF YOU USE IT
    double** _grad_g;   //NEED TO ALLOCATE THIS ONE BEFORE IF YOU USE IT
    double** _grad_g3D;  //for cross products   //NEED TO ALLOCATE THIS ONE BEFORE IF YOU USE IT
    double* _curl_g3D;   //NEED TO ALLOCATE THIS ONE BEFORE IF YOU USE IT
    std::vector < std::vector<double> > _el_average;/*[VB][spacedim]*/ //NEED TO ALLOCATE THIS EXPLICITLY WHERE IT'S USED... TODO this class must be reconsidered!!! with std::vectorss, and so on!!!
    uint _FEord; 
    uint _dim;
    uint _ndof[VB];
    
    EqnBase*  _eqnptr;
    Quantity* _qtyptr;
    
    CurrGaussPointBase & _currGP;
    CurrElem & _currEl;
    

    //TODO all these function are of the SET type (this is how I should call them), that is why they are NOT CONST
   void  VectWithQtyFillBasic();             //this needs the quantity and the fe map
   void                 val_g(const uint vbflag); //this only needs the CUrrent GAUSS  //No Quantity needed
   void                grad_g(const uint vbflag); //this only needs the CUrrent GAUSS  //No Quantity needed
   void                curl_g(const uint vbflag); //this only needs the CUrrent GAUSS  //No Quantity needed
   void            ExtendDofs(const uint vbflag); //this only needs the CUrrent GAUSS  //No Quantity needed
   void         GetElDofsVect(const uint vbflag, const uint Level); //this only needs the CUrrent ELEMENT
   void        SetElemAverage(const uint vbflag);
  
  //if you have NO Quantity and NO Equation ==========
//   void   SetElDofsFromArgs(const uint vb,const double * dofs);   //if you have NO Quantity and NO Equation, we should do the more flexible version of a Vect    
                                                                 //but the point is that we have to pass also the offset...
    
  };
  

#endif