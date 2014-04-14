/*=========================================================================

 Program: FEMUS
 Module: Parameter
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Parameter.hpp"


namespace femus {



Parameter::Parameter(const double Lref,const double Uref,const double DeltaTref) {
  _Lref = Lref;
  _Uref = Uref;
  _DeltaTref = DeltaTref;
}

double Parameter::Get_reference_lenght() {
  return _Lref;
}

double Parameter::Get_reference_velocity() {
  return _Uref;
}

double Parameter::Get_reference_temperature() {
  return _DeltaTref;
}




} //end namespace femus


