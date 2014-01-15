#include "Material.hpp"
#include "Parameter.hpp"

Material::Material(Parameter& par): _parameter(par) {
  _density = 1.;
  _thermal_conductivity = 1.;
  _heat_capacity = 1.;
  _thermal_expansion_coefficient = 1.e-06;
}

Material::Material(Parameter& par,const double density,const double k,
                   const double cp, const double alpha):
  _parameter(par) {
  _density = density;
  _thermal_conductivity = k;
  _heat_capacity = cp;
  _thermal_expansion_coefficient = alpha;
}

void Material::set_density(const double density) {
  _density = density;
}

double Material::get_density() {
  return _density;
}

double Material::get_thermal_conductivity() {
  return _thermal_conductivity;
}

double Material::get_heat_capacity() {
  return _heat_capacity;
}

double Material::get_thermal_expansion_coefficient() {
  return _thermal_expansion_coefficient;
}