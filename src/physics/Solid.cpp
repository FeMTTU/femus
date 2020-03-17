/*=========================================================================

 Program: FEMUS
 Module: Solid
 Authors: Simone Bnà, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <iostream>
#include <cstdlib>
#include <cstring>
#include "Solid.hpp"
#include "Material.hpp"
#include "Parameter.hpp"


namespace femus {



using namespace std;

Solid::Solid() : Material() {
  _young_module = 1.;
  _poisson_coeff = 0.3;
  _lambda_lame = 1.;
  _mu_lame  = 1.;
  _model = 0;
}

Solid::Solid(Parameter& par) : Material(par) {
  _young_module = 1.;
  _poisson_coeff = 0.3;
  _lambda_lame = 1.;
  _mu_lame  = 1.;
  _model = 0;
}

Solid::Solid(Parameter& par, 
             const double young_module,
             const double poisson_coeff,
             const double density, 
             const char model[],
             const double k,
             const double cp,
             const double alpha) : Material(par, density, k, cp, alpha) {
                 
  _young_module = young_module;


  if (!strcmp(model,"Linear_elastic") || !strcmp(model,"Saint-Venant")) {
    _model = 0;
    _penalty = false;
    _mass_penalty = false;
  }
  else if (!strcmp(model,"Saint-Venant-Penalty")) {
    _model = 0;
    _penalty = true;
    _mass_penalty = false;
  }
  else if (!strcmp(model,"Neo-Hookean") || !strcmp(model,"Neo-Hookean-MassPenalty")) {
    _model = 1;
    _penalty = false;
    _mass_penalty = !strcmp(model,"Neo-Hookean-MassPenalty") ? true : false;
  }
  else if (!strcmp(model,"Neo-Hookean-BW") || !strcmp(model,"Neo-Hookean-BW-MassPenalty") ) {  //Bonet-Wood
    _model = 2;
    _penalty = false;
    _mass_penalty = !strcmp(model,"Neo-Hookean-BW-MassPenalty") ? true : false;
  }
  else if (!strcmp(model,"Neo-Hookean-BW-Penalty")) {
    _model = 3;
    _penalty = true;
    _mass_penalty = false;
  }
  else if (!strcmp(model,"Neo-Hookean-AB-Penalty")) {  //Allan-Bower
    _model = 4;
    _penalty = true;
    _mass_penalty = false;
  }
   else if (!strcmp(model,"Mooney-Rivlin") || !strcmp(model,"Mooney-Rivlin-MassPenalty")) {
    _model = 5;
    _penalty = false;
    _mass_penalty = (!strcmp(model,"Mooney-Rivlin-MassPenalty")) ? true : false;
  }
  else {
    cout << "Error! This solid model is not implemented " << endl;
    abort();
  }

  if (poisson_coeff <= 0.5 && poisson_coeff >= 0) {
    _poisson_coeff = poisson_coeff;
  } 
  else {
    cout << "Error: the value for the Poisson coeffcient must be greater than 0 and less equal than 0.5!" << endl;
    abort();
  }
  
  
  if(poisson_coeff<0.5){
    _lambda_lame = (_young_module*_poisson_coeff)/((1.+_poisson_coeff)*(1.-2.*_poisson_coeff));
  }
  else if (true == _penalty){
    std::cout << "Error this solid model requires a Poisson coeffcient strictly less than 0.5"<<endl;
    abort();
  }
  else{
    cout << "Warning: the value for the Poisson coeffcient is 0.5, the material is incompressible"<<endl
	 << "The Lame constant is infinity and it has been set equal to 1.0e100" << endl;
    _lambda_lame = 1.0e100;
  }
  _mu_lame     = _young_module/(2.*(1.+_poisson_coeff));

//   cout << endl << "SOLID properties: " << endl;
//   cout << "Density [kg/m3]: " << _density << endl;
//   cout << "Young Module [Pa/m2]: " << _young_module << endl;
//   cout << "Poisson Coefficient: " << _poisson_coeff << endl;
//   cout << "Lambda (Lame first parameter): " << _lambda_lame << endl;
//   cout << "Shear Modulus (Lame second parameter): " << _mu_lame <<  endl << endl;

}

void Solid::set_young_module(const double young_module) {
  _young_module = young_module;
}

void Solid::set_poisson_coeff(const double poisson_coeff) {
  _poisson_coeff = poisson_coeff;
}

const unsigned Solid::get_physical_model() const{
  return _model;
}

const double Solid::get_young_module() const{
  return _young_module;
}

const double Solid::get_poisson_coeff() const{
  return _poisson_coeff;
}

const double Solid::get_lame_lambda() const{
  return _lambda_lame;
}

const double Solid::get_lame_shear_modulus() const{
  return _mu_lame;
}

const bool Solid::get_if_penalty() const{ 
  return _penalty;
}

const bool Solid::get_if_mass_penalty() const{ 
  return _mass_penalty;
}

std::ostream & operator << (std::ostream & os, const Solid & solid)
{
    
  os << "Density: [kg/m^3] "         << solid._density       << std::endl;
  os << "Young Module: [Pa] "        << solid._young_module  << std::endl;
  os << "Poisson coefficient: [Pa] " << solid._poisson_coeff << std::endl;
  os << "Lambda Lamé: [Pa] "         << solid._lambda_lame   << std::endl;
  os << "Mu lamé: [Pa] "             << solid._mu_lame       << std::endl;
  os << "Physical Model: "           << solid._model         << std::endl;
  os << std::endl;
  
  return os;
  
}

// Take a const-reference to the right-hand side of the assignment.
// Return a non-const reference to the left-hand side.
Solid& Solid::operator=(const Solid &solid) {
    
  this->_parameter 			= solid._parameter;
  this->_density 			= solid._density;
  this->_thermal_conductivity 		= solid._thermal_conductivity;
  this->_heat_capacity 			= solid._heat_capacity;
  this->_thermal_expansion_coefficient 	= solid._thermal_expansion_coefficient;
  this->_young_module 			= solid._young_module;
  this->_poisson_coeff 			= solid._poisson_coeff;
  this->_lambda_lame 			= solid._lambda_lame;
  this->_mu_lame 			= solid._mu_lame;
  this->_model 				= solid._model;
  this->_penalty 			= solid._penalty;
  this->_mass_penalty 			= solid._mass_penalty;
  
  return *this;  // Return a reference to myself.
  
}



} //end namespace femus


