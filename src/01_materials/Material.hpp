/*=========================================================================

 Program: FEMUS
 Module: Material
 Authors: Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_physics_Material_hpp__
#define __femus_physics_Material_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <iostream>


namespace femus {



//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class Parameter;

class Material {

public:

    /** costructor */
    Material(Parameter& par, const double density, const double k=1.,
             const double cp=1., const double alpha=1.e-06);

    /** simple costructor */
    Material(Parameter& par);

    /** simplest constructor */
    Material();

    /** destructor */
    ~Material() {};

    /** parameters for adimensionalization */
    Parameter* _parameter;

    /** set the density */
    void set_density(const double density);

    /** get the density */
    const double get_density() const;

    /** get the thermal conductivity */
    const double get_thermal_conductivity() const;

    /** get the heat capacity */
    const double get_heat_capacity() const;

    /** get the thermal expansion coefficient */
    const double get_thermal_expansion_coefficient() const;

    /** printing operator */
    friend std::ostream & operator << (std::ostream & os, const Material & mat);

    /** overloading of the = operator */
    Material & operator=(const Material &mat);

protected:

    double _density;

    double _thermal_conductivity;

    double _heat_capacity;

    double _thermal_expansion_coefficient;

};


} //end namespace femus



#endif
