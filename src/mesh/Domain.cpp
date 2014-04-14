#include "Domain.hpp"


namespace femus {





Domain::Domain(const uint spacedim_in, RunTimeMap<double> & map_in) :
  _spacedim(spacedim_in),_domain_rtmap(map_in) {  }



Domain::~Domain() {  }


} //end namespace femus


