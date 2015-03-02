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

#ifndef __femus_physics_Fluid_hpp__
#define __femus_physics_Fluid_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Material.hpp"


namespace femus {



//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class Parameter;

class Fluid : public Material {

public:

    /** constructors */
    Fluid(Parameter& par, const double viscosity, const double density, const char model[]="Newtonian",
          const double k=1., const double cp=1., const double alpha=1.e-06);

    Fluid(Parameter& par);

    Fluid();

    /** destructor */
    ~Fluid() {};

    /** To be Added */
    void set_viscosity(const double viscosity);

    /** To be Added */
    const double get_viscosity() const;

    /** To be Added */
    const double get_IReynolds_number() const;

    /** To be Added */
    const double get_Prandtl_number() const;

    /** To be Added */
    const double get_Rayleigh_number() const;

    /** To be Added */
    const unsigned get_physical_model() const;

    /** To be Added */
    const double get_Peclet_number() const;

    /** printing operator */
    friend std::ostream & operator << (std::ostream & os, const Fluid & fluid);

    /** overloading of the = operator */
    Fluid & operator=(const Fluid &fluid);

private:

    // member data
    double _viscosity; //< dynamic viscosity
    double _IRe;       //< Inverse Reynolds number
    double _Prandtl;   //< Prandtl number
    double _Froude;    //< Froude number
    double _Rayleigh;  //< Rayleigh number
    double _Peclet;    //< Peclet number
    double _Reynolds;  //< Reynolds number
    double _Grashof;   //< Grashof number
    unsigned _model;   //< Physical model for the fluid (0-->Newtonian)

};



} //end namespace femus



#endif
