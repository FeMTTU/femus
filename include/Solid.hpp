#ifndef __solid_hpp__
#define __solid_hpp__

#include "Material.hpp"

class Parameter;

class Solid : public Material {

private:
  double _young_module;
  double _poisson_coeff;
  double _lambda_lame;
  double _mu_lame;
  unsigned _model;

public:
  Solid(Parameter& par);
  Solid(Parameter& par, const double young_module, const double poisson_coeff,
        const double density, const char model[]= "Linear_elastic",
        const double k=1., const double cp=1., const double alpha=1.e-06);

  ~Solid() {};
  void set_young_module(const double young_module);
  void set_poisson_coeff(const double poisson_coeff);
  double get_young_module();
  double get_poisson_coeff();
  double get_lame_lambda();
  double get_lame_shear_modulus();
  unsigned get_physical_model();

};


#endif