/*=========================================================================

 Program: FEMUS
 Module: Material
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

#include "Material.hpp"
#include "Parameter.hpp"


namespace femus {



Material::Material() {
  _parameter = NULL;
  _density = 1.;
  _thermal_conductivity = 1.;
  _heat_capacity = 1.;
  _thermal_expansion_coefficient = 1.e-06;
}

Material::Material(Parameter& par) {
  _parameter = &par;
  _density = 1.;
  _thermal_conductivity = 1.;
  _heat_capacity = 1.;
  _thermal_expansion_coefficient = 1.e-06;
}

Material::Material(Parameter& par,const double density,const double k,
                   const double cp, const double alpha)
 {
  _parameter = &par;
  _density = density;
  _thermal_conductivity = k;
  _heat_capacity = cp;
  _thermal_expansion_coefficient = alpha;
}

void Material::set_density(const double density) {
  _density = density;
}

const double Material::get_density() const{
  return _density;
}

const double Material::get_thermal_conductivity() const{
  return _thermal_conductivity;
}

const double Material::get_heat_capacity() const{
  return _heat_capacity;
}

const double Material::get_thermal_expansion_coefficient() const{
  return _thermal_expansion_coefficient;
}

std::ostream & operator << (std::ostream & os, const Material & mat)
{
  os << "Density: " << mat.get_density() << std::endl;
  os << "Thermal conductivity: " << mat.get_thermal_conductivity() << std::endl;
  os << "Heat capacity: " << mat.get_heat_capacity() << std::endl;
  os << "Thermal expansion coefficient: " << mat.get_thermal_expansion_coefficient() << std::endl;
  os << std::endl;
  return os;
}

// Take a const-reference to the right-hand side of the assignment.
// Return a non-const reference to the left-hand side.
Material& Material::operator=(const Material &mat) {
  this->_parameter = mat._parameter;
  this->_density = mat._density;
  this->_thermal_conductivity = mat._thermal_conductivity;
  this->_heat_capacity = mat._heat_capacity;
  this->_thermal_expansion_coefficient = mat._thermal_expansion_coefficient;
  return *this;  // Return a reference to myself.
}



} //end namespace femus


