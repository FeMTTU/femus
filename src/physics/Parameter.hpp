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

#ifndef __parameter_hpp__
#define __parameter_hpp__


namespace femus {
  

class Parameter {

private:
  double _Lref;
  double _Uref;
  double _DeltaTref;

public:
  Parameter(const double Lref=1.,const double Uref=1.,const double DeltaTref=1.);
  double Get_reference_lenght();
  double Get_reference_velocity();
  double Get_reference_temperature();

};


} //end namespace femus


#endif