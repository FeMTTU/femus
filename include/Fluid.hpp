#ifndef __fluid_hpp__
#define __fluid_hpp__

#include "Material.hpp"

class Parameter;

class Fluid : public Material {

private:
  double _viscosity; // dynamic viscosity
  double _IRe;       // Inverse Reynolds number
  double _Prandtl;   // Prandtl number
  double _Froude;    // Froude number
  double _Rayleigh;  // Rayleigh number
  double _Peclet;    // Peclet number
  double _Reynolds;  // Reynolds number
  double _Grashof;   // Grashof number
  unsigned _model;   // Physical model for the fluid (0-->Newtonian)

public:
  Fluid(Parameter& par);
  Fluid(Parameter& par, const double viscosity, const double density, const char model[]="Newtonian",
        const double k=1., const double cp=1., const double alpha=1.e-06);
  ~Fluid() {};

  void set_viscosity(const double viscosity);
  double get_viscosity();
  double get_IReynolds_number();
  double get_Prandtl_number();
  double get_Rayleigh_number();
  unsigned get_physical_model();


};


#endif