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

#ifndef __fluid_hpp__
#define __fluid_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Material.hpp"

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
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
  
  Fluid(Parameter& par, const double viscosity, const double density, const char model[]="Newtonian",
        const double k=1., const double cp=1., const double alpha=1.e-06);
  
  Fluid(Parameter& par);
  
  Fluid();
  
  ~Fluid() {};

  void set_viscosity(const double viscosity);
  
  const double get_viscosity() const;
  
  const double get_IReynolds_number() const;
  
  const double get_Prandtl_number() const;
  
  const double get_Rayleigh_number() const;
  
  const unsigned get_physical_model() const;
  
  const double get_Peclet_number() const;
  
  /** printing operator */
  friend std::ostream & operator << (std::ostream & os, const Fluid & fluid); 

  /** overloadinf of the = operator */
  Fluid & operator=(const Fluid &fluid);

};


#endif