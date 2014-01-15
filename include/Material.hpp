#ifndef __material_hpp__
#define __material_hpp__

class Parameter;

class Material {

protected:
  double _density;
  double _thermal_conductivity;
  double _heat_capacity;
  double _thermal_expansion_coefficient;

public:
  Material(Parameter& par, const double density, const double k=1.,
           const double cp=1., const double alpha=1.e-06);
  Material(Parameter& par);
  ~Material() {};
  Parameter & _parameter;
  void set_density(const double density);
  double get_density();
  double get_thermal_conductivity();
  double get_heat_capacity();
  double get_thermal_expansion_coefficient();

};

#endif