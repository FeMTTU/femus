#include "Fluid.hpp"
#include "Material.hpp"
#include "iostream"
#include "cstdlib"
#include "cstring"
#include "Parameter.hpp"
using namespace std;


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

  const  double lref = par.Get_reference_lenght();
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

  cout << endl << "FLUID properties: " << endl;
  cout << "Density [kg/m3]: " << density << endl;
  cout << "Viscosity [Pa*s]: " << viscosity << endl;
  cout << "Heat capacity [J/Kg*K]: " << cp << endl;
  cout << "Thermal conductivity [W/m*K]: " << k << endl;
  cout << "Thermal expansion coefficient [1/K]: " << alpha << endl;
  cout << "Reynolds Number: " << _Reynolds << endl;
  cout << "Prandtl number: " << _Prandtl << endl;
  cout << "Froude Number: " << _Froude << endl;
  cout << "Grashof Number: " << _Grashof << endl;
  cout << "Rayleigh Number: " << _Rayleigh << endl;
  cout << "Peclet Number: " << _Peclet << endl << endl;
}

void Fluid::set_viscosity(const double viscosity) {
  _viscosity = viscosity;
}

double Fluid::get_viscosity() {
  return _viscosity;
}

unsigned Fluid::get_physical_model() {
  return _model;
}

double Fluid::get_IReynolds_number() {
  return _IRe;
}

double Fluid::get_Prandtl_number() {
  return _Prandtl;
};

double Fluid::get_Rayleigh_number() {
  return _Rayleigh;
};

double Fluid::get_Peclet_number() {
  return _Peclet;
};

