#include "Domain.hpp"

#include "Utils.hpp"




Domain::Domain(Utils & utils_in) :
 _utils(utils_in),_spacedim( (uint) _utils._urtmap.get("dimension")) {


}



Domain::~Domain() {  }