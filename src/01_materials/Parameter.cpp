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

double Parameter::Get_reference_length() {
  return _Lref;
}

double Parameter::Get_reference_velocity() {
  return _Uref;
}

double Parameter::Get_reference_temperature() {
  return _DeltaTref;
}


Gravity::Gravity(const double gx,const double gy,const double gz) {
  _g[0] = gx;
  _g[1] = gy;
  _g[2] = gz;
}

const double * Gravity::get_values() const {
  return _g;
}

std::ostream & operator << (std::ostream & os, const Gravity & g)
{
  os << "gx: " << g.get_values()[0] << std::endl;
  os << "gy: " << g.get_values()[1] << std::endl;
  os << "gz: " << g.get_values()[2] << std::endl;
  os << std::endl;
  return os;
}

// Take a const-reference to the right-hand side of the assignment.
// Return a non-const reference to the left-hand side.
Gravity& Gravity::operator=(const Gravity &g) {
  this->_g[0] = g._g[0];
  this->_g[1] = g._g[1];
  this->_g[2] = g._g[2];
  return *this;  // Return a reference to myself.
}


} //end namespace femus


