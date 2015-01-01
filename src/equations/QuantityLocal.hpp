#ifndef __quantitylocal_h__
#define __quantitylocal_h__

#include <vector>
#include "VBTypeEnum.hpp"
#include "CurrGaussPoint.hpp"


namespace femus {

 
 class EqnBase;
 class Quantity;
 class CurrGaussPointBase;
 class CurrElem;

  class QuantityLocal {
    
  public:
    
     QuantityLocal(CurrGaussPointBase &, CurrElem &);
    ~QuantityLocal();

    
    //TODO all these function are of the SET type (this is how I should call them), that is why they are NOT CONST
   void  VectWithQtyFillBasic();             //this needs the quantity and the fe map
   void Allocate();
   void Deallocate();
   void                 val_g(); //this only needs the CUrrent GAUSS  //No Quantity needed
   void                grad_g(); //this only needs the CUrrent GAUSS  //No Quantity needed
   void                curl_g(); //this only needs the CUrrent GAUSS  //No Quantity needed
   void            ExtendDofs(); //this only needs the CUrrent GAUSS  //No Quantity needed
   void         GetElDofsVect(const uint Level); //this only needs the CUrrent ELEMENT
   void        SetElemAverage();
  
  //if you have NO Quantity and NO Equation ==========
//   void   SetElDofsFromArgs(const uint vb,const double * dofs);   //if you have NO Quantity and NO Equation, we should do the more flexible version of a Vect    
                                                                 //but the point is that we have to pass also the offset...
    

    double*  _val_g;                                      //NEED TO ALLOCATE THIS ONE BEFORE IF YOU USE IT  //WHY ARE WE NOT ALLOCATING IN THE CONSTRUCTOR?!?
    double*  _val_g3D;   //for cross products             //NEED TO ALLOCATE THIS ONE BEFORE IF YOU USE IT
    double*  _val_dofs;   //NEED TO ALLOCATE THIS ONE BEFORE IF YOU USE IT
    double*  _val_dofs3D;   //NEED TO ALLOCATE THIS ONE BEFORE IF YOU USE IT
    double** _grad_g;   //NEED TO ALLOCATE THIS ONE BEFORE IF YOU USE IT
    double** _grad_g3D;  //for cross products   //NEED TO ALLOCATE THIS ONE BEFORE IF YOU USE IT
    double* _curl_g3D;   //NEED TO ALLOCATE THIS ONE BEFORE IF YOU USE IT
    
    std::vector<double> _el_average;  /*[spacedim]*/ //NEED TO ALLOCATE THIS EXPLICITLY WHERE IT'S USED... TODO this class must be reconsidered!!! with std::vectorss, and so on!!!

    uint _FEord; 
    uint _dim;
    uint _ndof;
    Quantity* _qtyptr;
    EqnBase*  _eqnptr;
    
    inline const CurrElem &  GetCurrentElem() const { 
      return _currEl;
    }
    
  protected:
        
    CurrGaussPointBase & _currGP;
    CurrElem & _currEl;
    
  };
  


} //end namespace femus



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
 
 //Ok, so Vect can have a quantity with equation,
 //or a quantity without equation,
 //or no quantity at all
 
#endif