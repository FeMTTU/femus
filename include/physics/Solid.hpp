/*=========================================================================

 Program: FEMUS
 Module: FemTTUInit
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __solid_hpp__
#define __solid_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Material.hpp"


namespace femus {



//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
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
  
  Solid();

  ~Solid() {};
  
  void set_young_module(const double young_module);
  
  void set_poisson_coeff(const double poisson_coeff);
  
  const double get_young_module() const;
  
  const double get_poisson_coeff() const;
  
  const double get_lame_lambda() const;
  
  const double get_lame_shear_modulus() const;
  
  const unsigned get_physical_model() const;
  
  /** printing operator */
  friend std::ostream & operator << (std::ostream & os, const Solid & solid); 

  /** overloadinf of the = operator */
  Solid & operator=(const Solid &solid);

};



} //end namespace femus



#endif