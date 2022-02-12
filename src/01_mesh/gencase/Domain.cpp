/*=========================================================================

 Program: FEMUS
 Module: Domain
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "Domain.hpp"


namespace femus {





Domain::Domain(const uint spacedim_in, const FemusInputParser<double> & map_in) :
  _spacedim(spacedim_in),_domain_rtmap(map_in) {  }



Domain::~Domain() {  }


} //end namespace femus


