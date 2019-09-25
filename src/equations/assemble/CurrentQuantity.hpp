#ifndef __femus_equations_CurrentQuantity_hpp__
#define __femus_equations_CurrentQuantity_hpp__

#include <vector>
#include "CurrentElem.hpp"



namespace femus {

 
//  class CurrentElem<double>; 
 class SystemTwo;
 class Quantity;

 //Remember that you need to allocate the operators before if you use them
 
  class CurrentQuantity {
    
  public:
    
     CurrentQuantity(const CurrentElem<double> & curr_el_in);

    
    //TODO all these function are of the SET type (this is how I should call them), that is why they are NOT CONST
   void  VectWithQtyFillBasic();             //this needs the quantity and the fe map
   void Allocate();

   void                 val_g(); //this only needs the CUrrent GAUSS  //No Quantity needed
   void                grad_g(); //this only needs the CUrrent GAUSS  //No Quantity needed
   void                curl_g(); //this only needs the CUrrent GAUSS  //No Quantity needed
   void            ExtendDofs(); //this only needs the CUrrent GAUSS  //No Quantity needed
   void         GetElemDofs(); //this only needs the CUrrent ELEMENT
//    void         GetElemDofs_two(); //this only needs the CUrrent ELEMENT
   void         GetElemDofs(const std::vector<CurrentQuantity*> vec_in);
   void        SetElemAverage();
  
  //if you have NO Quantity and NO Equation ==========
//   void   SetElDofsFromArgs(const uint vb,const double * dofs);   //if you have NO Quantity and NO Equation, we should do the more flexible version of a Vect    
                                                                 //but the point is that we have to pass also the offset...

    std::vector<double>  _val_g;
    std::vector<double>  _val_g3D;
    std::vector< std::vector<double> >  _grad_g;
    std::vector< std::vector<double> >  _grad_g3D;
    std::vector<double> _curl_g3D;
    
    std::vector<double>  _val_dofs;  
    std::vector<double>  _val_dofs3D;
    
    std::vector<double> _el_average;
    
    Quantity* _qtyptr;
    SystemTwo*  _eqnptr;
    
    inline const CurrentElem<double> &  GetCurrentElem() const; 
    
    
    std::string _SolName;
    
    vector< int > _dofSOL;
    unsigned _indexSOL;
    unsigned _indSOL;
    unsigned _SolTypeSOL;

// once you set these three, you can call "Allocate"    
    uint _FEord; 
    uint _dim;
    uint _ndof;
// once you set these three, you can call "Allocate"    
    
  protected:
        
    const CurrentElem<double> & _currEl;
    
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
