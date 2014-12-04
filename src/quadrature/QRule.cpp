#include "QRule.hpp"

#include <cstdlib>
#include <iostream>
#include <cstring>

#include "GeomElTypeEnum.hpp"
#include "GeomEl.hpp"


namespace femus {

 
 QRule::QRule(const GeomEl & geomel_in)  : _geomel(geomel_in),_qrule_type("fifth")  {

   
      if (!strcmp(_geomel._geomel_id.c_str(),"hex")) {  
		      _NoGaussVB = 27; 
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"wedge")) { 
	std::cout << "No gauss points on wedge" << std::endl; abort();
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"tet")) {
		      _NoGaussVB = 15; 
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"quad")) {
      		      _NoGaussVB = 9;
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"tri")) {
      		      _NoGaussVB = 7;
      }	
      else if (!strcmp(_geomel._geomel_id.c_str(),"line")) {
		      _NoGaussVB = 3;
      }
      else {
	std::cout << _geomel._geomel_id << " is not a valid option" << std::endl; 
	abort();
      }
       
//now fill the weights =========================
    _weightVB.resize(_NoGaussVB);
       

      if (!strcmp(_geomel._geomel_id.c_str(),"hex")) {  

_weightVB[0] =    0.17146776406036;
_weightVB[1] =    0.27434842249657;  
_weightVB[2] =    0.17146776406036;  
_weightVB[3] =    0.27434842249657;
_weightVB[4] =    0.43895747599451;
_weightVB[5] =    0.27434842249657;
_weightVB[6]  =   0.17146776406036;
_weightVB[7]  =   0.27434842249657;
_weightVB[8]  =   0.17146776406036;
_weightVB[9]  =   0.27434842249657;
_weightVB[10] =   0.43895747599451;
_weightVB[11] =   0.27434842249657;
_weightVB[12] =   0.43895747599451;
_weightVB[13] =   0.70233196159122;
_weightVB[14] =   0.43895747599451;
_weightVB[15] =   0.27434842249657;
_weightVB[16] =   0.43895747599451;
_weightVB[17] =   0.27434842249657;
_weightVB[18] =   0.17146776406036;
_weightVB[19] =   0.27434842249657;
_weightVB[20] =   0.17146776406036;
_weightVB[21] =   0.27434842249657;
_weightVB[22] =   0.43895747599451;
_weightVB[23] =   0.27434842249657;
_weightVB[24] =   0.17146776406036;
_weightVB[25] =   0.27434842249657;
_weightVB[26] =   0.17146776406036;

      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"wedge")) { 
	std::cout << "No gauss points on wedge" << std::endl; abort();
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"tet")) {
_weightVB[0] =    0.030283678097089;
_weightVB[1] =    0.006026785714286;  
_weightVB[2] =    0.006026785714286;  
_weightVB[3] =    0.006026785714286;
_weightVB[4] =    0.006026785714286;
_weightVB[5] =    0.011645249086029;
_weightVB[6]  =   0.011645249086029;
_weightVB[7]  =   0.011645249086029;
_weightVB[8]  =   0.011645249086029;
_weightVB[9]  =   0.010949141561386;
_weightVB[10] =   0.010949141561386;
_weightVB[11] =   0.010949141561386;
_weightVB[12] =   0.010949141561386;
_weightVB[13] =   0.010949141561386;     
_weightVB[14] =   0.010949141561386;
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"quad")) {

_weightVB[0] =   0.30864197530864    ;
_weightVB[1] =   0.49382716049383    ;
_weightVB[2] =   0.30864197530864    ;
_weightVB[3] =   0.49382716049383    ;
_weightVB[4] =   0.79012345679012    ;
_weightVB[5] =   0.49382716049383    ;
_weightVB[6]  =  0.30864197530864    ;
_weightVB[7]  =  0.49382716049383    ;
_weightVB[8]  =  0.30864197530864    ;

      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"tri")) {
_weightVB[0] =   0.1125;  
_weightVB[1] =   0.062969590272414;  
_weightVB[2] =   0.062969590272414;  
_weightVB[3] =   0.062969590272414;  
_weightVB[4] =   0.066197076394253;  
_weightVB[5] =   0.066197076394253;  
_weightVB[6]  =  0.066197076394253;  
      }	
      else if (!strcmp(_geomel._geomel_id.c_str(),"line")) {
_weightVB[0] =   0.555555555555555555555555555556    ;
_weightVB[1] =   0.888888888888888888888888888889    ;
_weightVB[2] =   0.555555555555555555555555555556    ;
      }
      else {
	std::cout << _geomel._geomel_id << " is not a valid option" << std::endl; 
	abort();
      }
  
  
       
 }

 


} //end namespace femus
