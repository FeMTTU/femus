/*=========================================================================

 Program: FEMUS
 Module: Parameter
 Authors: Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_physics_Parameter_hpp__
#define __femus_physics_Parameter_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <iostream>


namespace femus {


class Parameter {

public:

    /** constructor */
    Parameter(const double Lref=1.,const double Uref=1.,const double DeltaTref=1.);

    /** To be Added */
    double Get_reference_length();

    /** To be Added */
    double Get_reference_velocity();

    /** To be Added */
    double Get_reference_temperature();

private:

    // member data

    double _Lref;
    double _Uref;
    double _DeltaTref;

};

class Gravity {

public:

    /** constructor */
    Gravity(const double gx=0.,const double gy=0.,const double gz=0.);

    /** get an array containing gx, gy and gx */
    const double * get_values() const;
    
    /** printing operator */
    friend std::ostream & operator << (std::ostream & os, const Gravity & g);

    /** overloading of the = operator */
    Gravity & operator=(const Gravity &g);

private:

    // member data

    double _g[3];

};


} //end namespace femus


#endif
