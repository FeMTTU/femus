/*=========================================================================

 Program: FEMUS
 Module: Fluid
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
#include "Fluid.hpp"
#include "Material.hpp"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include "Parameter.hpp"


namespace femus {


using namespace std;

Fluid::Fluid() : Material() {
  _viscosity = 1.;
  _IRe = 1;
  _Reynolds = 1.;
  _model = 0;
  _Prandtl=1.;
  _Froude = 1.;
  _Rayleigh = 1.;
  _Peclet = 1.;
  _Grashof = 1.;
}

Fluid::Fluid(Parameter& par) : Material(par) {
  _viscosity = 1.;
  _IRe = 1;
  _Reynolds = 1.;
  _model = 0;
  _Prandtl=1.;
  _Froude = 1.;
  _Rayleigh = 1.;
  _Peclet = 1.;
  _Grashof = 1.;
}

Fluid::Fluid(Parameter& par, const double viscosity, const double density, const char model[],
             const double k, const double cp, const double alpha): Material(par,density,k,cp, alpha) {
//   cout << "calling 2" << endl;

  const  double lref = par.Get_reference_length();
  const double uref = par.Get_reference_velocity();
  const double DeltaTref = par.Get_reference_temperature();
  _viscosity = viscosity;
  _density   = density;
  _Reynolds = (_density*uref*lref)/_viscosity;
  _IRe = 1./_Reynolds;  // rho*U*Lref/mu

  if (!strcmp(model,"Newtonian")) {
    _model = 0;
  } else {
    cout<<"Error! This fluid model is not implemented "<<endl;
    exit(1);
  }

  _Prandtl = (cp*viscosity)/k;
  _Froude = (uref*uref)/(9.80665*lref);    // U^2/Lref*g
  _Grashof = (density*density*9.80665*alpha*DeltaTref*lref*lref*lref)/(viscosity*viscosity);  //*L^3 * delta_T

  _Rayleigh = _Grashof*_Prandtl;
  _Peclet = _Prandtl*_Reynolds;

//   cout << endl << "FLUID properties: " << endl;
//   cout << "Density [kg/m3]: " << density << endl;
//   cout << "Viscosity [Pa*s]: " << viscosity << endl;
//   cout << "Heat capacity [J/Kg*K]: " << cp << endl;
//   cout << "Thermal conductivity [W/m*K]: " << k << endl;
//   cout << "Thermal expansion coefficient [1/K]: " << alpha << endl;
//   cout << "Reynolds Number: " << _Reynolds << endl;
//   cout << "Prandtl number: " << _Prandtl << endl;
//   cout << "Froude Number: " << _Froude << endl;
//   cout << "Grashof Number: " << _Grashof << endl;
//   cout << "Rayleigh Number: " << _Rayleigh << endl;
//   cout << "Peclet Number: " << _Peclet << endl << endl;
}

void Fluid::set_viscosity(const double viscosity) {
  _viscosity = viscosity;
}

const double Fluid::get_viscosity() const {
  return _viscosity;
}

const unsigned Fluid::get_physical_model() const{
  return _model;
}

const double Fluid::get_IReynolds_number() const{
  return _IRe;
}

const double Fluid::get_Prandtl_number() const{
  return _Prandtl;
};

const double Fluid::get_Rayleigh_number() const{
  return _Rayleigh;
};

const double Fluid::get_Peclet_number() const{
  return _Peclet;
};

std::ostream & operator << (std::ostream & os, const Fluid & fluid)
{
  os << "Density: [kg/m3] " << fluid._density << std::endl;
  os << "viscosity: [Pa*s] " << fluid._viscosity << std::endl;
  os << "Inverse of Reynolds number: " << fluid._IRe << std::endl;
  os << "Prandtl number: " << fluid._Prandtl << std::endl;
  os << "Froude number: " << fluid._Froude << std::endl;
  os << "Rayleigh number: " << fluid._Rayleigh << std::endl;
  os << "Peclet number: " << fluid._Peclet << std::endl;
  os << "Reynolds number: " << fluid._Reynolds << std::endl;
  os << "Grashof number: " << fluid._Grashof << std::endl;
  os << "Physical Model: " << fluid._model << std::endl;
  os << std::endl;
  
  return os;
}

// Take a const-reference to the right-hand side of the assignment.
// Return a non-const reference to the left-hand side.
Fluid& Fluid::operator=(const Fluid &fluid) {
    
  this->_parameter = fluid._parameter;
  this->_density = fluid._density;
  this->_thermal_conductivity = fluid._thermal_conductivity;
  this->_heat_capacity = fluid._heat_capacity;
  this->_thermal_expansion_coefficient = fluid._thermal_expansion_coefficient;
  this->_viscosity = fluid._viscosity;
  this->_IRe = fluid._IRe;
  this->_Prandtl = fluid._Prandtl;
  this->_Froude = fluid._Froude;
  this->_Rayleigh = fluid._Rayleigh;
  this->_Peclet = fluid._Peclet;
  this->_Reynolds = fluid._Reynolds;
  this->_Grashof = fluid._Grashof;
  this->_model = fluid._model;
  
  return *this;  // Return a reference to myself.
  
}



} //end namespace femus


