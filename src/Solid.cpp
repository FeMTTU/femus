#include "Solid.hpp"
#include "Material.hpp"
#include "iostream"
#include "cstdlib"
#include "cstring"
#include "Parameter.hpp"

using namespace std;

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

  cout << endl << "SOLID properties: " << endl;
  cout << "Density [kg/m3]: " << _density << endl;
  cout << "Young Module [Pa/m2]: " << _young_module << endl;
  cout << "Poisson Coefficient: " << _poisson_coeff << endl;
  cout << "Lambda (Lame first parameter): " << _lambda_lame << endl;
  cout << "Shear Modulus (Lame second parameter): " << _mu_lame <<  endl << endl;

}

unsigned Solid::get_physical_model() {
  return _model;
}


void Solid::set_young_module(const double young_module) {
  _young_module = young_module;
}

void Solid::set_poisson_coeff(const double poisson_coeff) {
  _poisson_coeff = poisson_coeff;
}

double Solid::get_young_module() {
  return _young_module;
}

double Solid::get_poisson_coeff() {
  return _poisson_coeff;
}

double Solid::get_lame_lambda() {
  return _lambda_lame;
}

double Solid::get_lame_shear_modulus() {
  return _mu_lame;
}
