/*=========================================================================

 Program: FEMUS
 Module: Quantity
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "Quantity.hpp"

#include "MultiLevelMeshTwo.hpp"

#include <vector>


namespace femus {



//What does the Quantity need? 
//The Physics? Yes, for the ref values, and so on
//The Mesh? Yes,through the physics
//The FE values?No
//The Equation pointers? Yes. For every quantity i'll have to say
//what is the equation that computes them, if there is one.
//I'll pass the equation pointer not directly to the constructor,
//but after I have the corresponding equation: set_eqn() or sthg

//If there is an equation, i'll put it; otherwise, i'll put NULL
//If this physical quantity is not considered at all for my simulation,
//then I won't have to instantiate it at all
//In the same way, if an equation is not considered, then you dont even
//instantiate it
//that is why you put things within TEMP_QTY and so on
//Clearly, the quantities that are relevant are dependent 
//of the EQUATIONS
//So, since you have b_hom and Bext, you decide to define two separate
//quantities MagFieldHom and MagFieldExt, because each of the two is 
//naturally attached to a different equation.
//For Pressure, for instance, one may have either NS equation,
//or pressure single projection equation.
//What you only have to change is the GetElDofs() function:
//if you pick from the pressure equation you pick the dofs at some point,
//if you pick from the NS equation you pick the dofs in another point.
//Basically what changes is the DofMap! So, taken one equation,
//the GetElDof function will have to be able of picking 
//either the QUADRATIC or the LINEAR variable.
//For now, we cannot pick only one component, but we can pick
//either ALL the QUADRATIC, or ALL the LINEAR components
//so, the GetElDofs is a function belonging to a specific EQUATION;
//so, we can make it VIRTUAL, so that every Eqn will have its own way
//of implementing it.

  Quantity::Quantity(std::string name_in,QuantityMap& qtymap_in, uint dim_in, uint FEord_in) : _qtymap(qtymap_in) {

   _name  = name_in;
   _dim   = dim_in;
   _FEord = FEord_in;
   _eqn   = NULL; //initialize to NULL: by default the Quantity has no associated Equation
   _refvalue.resize(_dim);
   for (uint i=0; i<_dim; i++)  _refvalue[i] = 1.; //default
   
   
  }

 Quantity::~Quantity() {}



  void Quantity::set_eqn(SystemTwo* eqn_in)  { 

    _eqn = eqn_in;
    
    return; 
  }

} //end namespace femus


