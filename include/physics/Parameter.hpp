#ifndef __parameter_hpp__
#define __parameter_hpp__

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

#endif