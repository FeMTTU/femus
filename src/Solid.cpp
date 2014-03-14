#include "Solid.hpp"
#include "Material.hpp"
#include "iostream"
#include "cstdlib"
#include "cstring"
#include "Parameter.hpp"

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

Solid::Solid(Parameter& par, const double young_module, const double poisson_coeff,
             const double density, const char model[],
             const double k, const double cp, const double alpha) : Material(par,density,k,cp,alpha) {
  _young_module = young_module;


  if (!strcmp(model,"Linear_elastic")) {
    _model = 0;
  } else if (!strcmp(model,"Neo-Hookean")) {
    _model = 1;
  } else {
    cout<<"Error! This solid model is not implemented "<<endl;
    exit(1);
  }

  if (poisson_coeff < 0.5 && poisson_coeff >= 0) {
    _poisson_coeff = poisson_coeff;
  } else {
    cout << "Error: the value for the Poisson coeffcient must be greater than 0 and lower than 0.5!" << endl;
    exit(1);
  }

  _lambda_lame = (_young_module*_poisson_coeff)/((1.+_poisson_coeff)*(1.-2.*_poisson_coeff));
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

std::ostream & operator << (std::ostream & os, const Solid & solid)
{
  os << "Density: " << solid._density << std::endl;
  os << "Young Module: " << solid._young_module << std::endl;
  os << "Poisson coeffcient: " << solid._poisson_coeff << std::endl;
  os << "Lambda Lamé: " << solid._lambda_lame << std::endl;
  os << "Mu lamé: " << solid._mu_lame << std::endl;
  os << "Physical Model: " << solid._model << std::endl;
  os << std::endl;
  return os;
}

// Take a const-reference to the right-hand side of the assignment.
// Return a non-const reference to the left-hand side.
Solid& Solid::operator=(const Solid &solid) {
  this->_parameter = solid._parameter;
  this->_density = solid._density;
  this->_thermal_conductivity = solid._thermal_conductivity;
  this->_heat_capacity = solid._heat_capacity;
  this->_thermal_expansion_coefficient = solid._thermal_expansion_coefficient;
  this->_young_module = solid._young_module;
  this->_poisson_coeff = solid._poisson_coeff;
  this->_lambda_lame = solid._lambda_lame;
  this->_mu_lame = solid._mu_lame;
  this->_model = solid._model;
  return *this;  // Return a reference to myself.
}