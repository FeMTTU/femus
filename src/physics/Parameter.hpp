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


namespace femus {


class Parameter {

public:

    /** constructor */
    Parameter(const double Lref=1.,const double Uref=1.,const double DeltaTref=1.);

    /** To be Added */
    double Get_reference_lenght();

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


} //end namespace femus


#endif
