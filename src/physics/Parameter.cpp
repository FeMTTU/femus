#include "Parameter.hpp"

Parameter::Parameter(const double Lref,const double Uref,const double DeltaTref) {
  _Lref = Lref;
  _Uref = Uref;
  _DeltaTref = DeltaTref;
}

double Parameter::Get_reference_lenght() {
  return _Lref;
}

double Parameter::Get_reference_velocity() {
  return _Uref;
}

double Parameter::Get_reference_temperature() {
  return _DeltaTref;
}

