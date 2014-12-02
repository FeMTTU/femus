#include "QRule.hpp"

#include <cstdlib>
#include <iostream>
 
#include "GeomElTypeEnum.hpp"
#include "GeomEl.hpp"


namespace femus {


 

  QRule::~QRule()  {
    
    delete  []  _weightVB;
 
 }
 
 
 QRule::QRule(GeomEl geomel_in)  : _geomel(geomel_in),_qrule_type("Gauss5th")  {
   
   
       switch(_geomel._geomel_type) {
	 case(QUADR) : {

               switch(_geomel._dim) {
		    
		    case(3) : {
		      _NoGaussVB = 27; 
		      break;
		    }

		    case(2) : {
      		      _NoGaussVB = 9;
		      break;
		    }

		    case(1) : {
		      _NoGaussVB = 3;
		      break;
		    }
		      
		    default: {std::cout << "FE: Dimension not implemented" << std::endl; abort();}
		    
		  } //end dim
	   
	   break;
         } //end case quadrilateral
      
	 case(TRIANG) :  {
               switch(_geomel._dim) {
		    
		    case(3) : {
		      _NoGaussVB = 15; 
		      break;
		    }

		    case(2) : {
      		      _NoGaussVB = 7;
		      break;
		    }

		    case(1) : {
		      _NoGaussVB = 3;
		      break;
		    }
		    
		    default: {std::cout << "FE: Dimension not implemented" << std::endl; abort();}
		    
		  } //end dim 

	   break;
          } //end case triangular
      
       } //end switch geomel
       
       
       
//now fill the weights =========================

    _weightVB = new double[_NoGaussVB];
       


        switch(_geomel._geomel_type) {
	 case(QUADR) : {

               switch(_geomel._dim) {
		    
		    case(3) : {

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
		      
		      break;
		    }

		    case(2) : {

_weightVB[0] =   0.30864197530864    ;
_weightVB[1] =   0.49382716049383    ;
_weightVB[2] =   0.30864197530864    ;
_weightVB[3] =   0.49382716049383    ;
_weightVB[4] =   0.79012345679012    ;
_weightVB[5] =   0.49382716049383    ;
_weightVB[6]  =  0.30864197530864    ;
_weightVB[7]  =  0.49382716049383    ;
_weightVB[8]  =  0.30864197530864    ;

		      break;
		    }
		    
		    case(1) : {
		      
_weightVB[0] =   0.555555555555555555555555555556    ;
_weightVB[1] =   0.888888888888888888888888888889    ;
_weightVB[2] =   0.555555555555555555555555555556    ;

                      break;
		    }

		    default: {std::cout << "FE: Dimension not implemented" << std::endl; abort();}
		    
		  } //end dim
	   
	   break;
         } //end case quadrilateral
      
	 case(TRIANG) :  {
               switch(_geomel._dim) {
		    
		    case(3) : {
		      
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
		      
		      
		      break;
		    }

		    case(2) : {
		      
_weightVB[0] =   0.1125;  
_weightVB[1] =   0.062969590272414;  
_weightVB[2] =   0.062969590272414;  
_weightVB[3] =   0.062969590272414;  
_weightVB[4] =   0.066197076394253;  
_weightVB[5] =   0.066197076394253;  
_weightVB[6]  =  0.066197076394253;  

		      break;
		    }
		    
    		    case(1) : {
		      
_weightVB[0] =   0.555555555555555555555555555556    ;
_weightVB[1] =   0.888888888888888888888888888889    ;
_weightVB[2] =   0.555555555555555555555555555556    ;

		      break;
		    }


		    default: {std::cout << "FE: Dimension not implemented" << std::endl; abort();}
		    
		  } //end dim 

	   break;
          } //end case triangular
      
       } //end switch geomel   
  
  
       
 }
 
 
 


} //end namespace femus


 