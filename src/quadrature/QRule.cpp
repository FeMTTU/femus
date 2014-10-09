#include "QRule.hpp"

#include <cstdlib>
#include <iostream>
 
#include "GeomElTypeEnum.hpp"
#include "GeomEl.hpp"


namespace femus {


 

  QRule::~QRule()  {
    
 for (uint vb=0; vb < VB; vb++)  {
    delete []        _weightVB[vb];
  }
  
 }
 
 
 QRule::QRule(GeomEl* geomel_in)  : _geomel(geomel_in),_qrule_type("Gauss5th")  {
   
   
       switch(_geomel->_geomel_type) {
	 case(QUADR) : {

               switch(_geomel->_dim) {
		    
		    case(3) : {
		      _NoGaussVB[VV] = 27;
		      _NoGaussVB[BB] = 9; 
		      break;
		    }

		    case(2) : {
      		      _NoGaussVB[VV] = 9;
		      _NoGaussVB[BB] = 3;
		      break;
		    }

		    default: {std::cout << "FE: Dimension not implemented" << std::endl; abort();}
		    
		  } //end dim
	   
	   break;
         } //end case quadrilateral
      
	 case(TRIANG) :  {
               switch(_geomel->_dim) {
		    
		    case(3) : {
		      _NoGaussVB[VV] = 15;
		      _NoGaussVB[BB] = 7; 
		      break;
		    }

		    case(2) : {
      		      _NoGaussVB[VV] = 7;
		      _NoGaussVB[BB] = 3;
		      break;
		    }

		    default: {std::cout << "FE: Dimension not implemented" << std::endl; abort();}
		    
		  } //end dim 

	   break;
          } //end case triangular
      
       } //end switch geomel
       
       
       
//now fill the weights =========================

for (uint vb=0; vb < VB; vb++)    _weightVB[vb] = new double[_NoGaussVB[vb]];
       


        switch(_geomel->_geomel_type) {
	 case(QUADR) : {

               switch(_geomel->_dim) {
		    
		    case(3) : {

// ========= VV ========
                                           
_weightVB[VV][0] =    0.17146776406036;
_weightVB[VV][1] =    0.27434842249657;  
_weightVB[VV][2] =    0.17146776406036;  
_weightVB[VV][3] =    0.27434842249657;
_weightVB[VV][4] =    0.43895747599451;
_weightVB[VV][5] =    0.27434842249657;
_weightVB[VV][6]  =   0.17146776406036;
_weightVB[VV][7]  =   0.27434842249657;
_weightVB[VV][8]  =   0.17146776406036;
_weightVB[VV][9]  =   0.27434842249657;
_weightVB[VV][10] =   0.43895747599451;
_weightVB[VV][11] =   0.27434842249657;
_weightVB[VV][12] =   0.43895747599451;
_weightVB[VV][13] =   0.70233196159122;
_weightVB[VV][14] =   0.43895747599451;
_weightVB[VV][15] =   0.27434842249657;
_weightVB[VV][16] =   0.43895747599451;
_weightVB[VV][17] =   0.27434842249657;
_weightVB[VV][18] =   0.17146776406036;
_weightVB[VV][19] =   0.27434842249657;
_weightVB[VV][20] =   0.17146776406036;
_weightVB[VV][21] =   0.27434842249657;
_weightVB[VV][22] =   0.43895747599451;
_weightVB[VV][23] =   0.27434842249657;
_weightVB[VV][24] =   0.17146776406036;
_weightVB[VV][25] =   0.27434842249657;
_weightVB[VV][26] =   0.17146776406036;
		      
// ========= BB ========
_weightVB[BB][0] =   0.30864197530864    ;
_weightVB[BB][1] =   0.49382716049383    ;
_weightVB[BB][2] =   0.30864197530864    ;
_weightVB[BB][3] =   0.49382716049383    ;
_weightVB[BB][4] =   0.79012345679012    ;
_weightVB[BB][5] =   0.49382716049383    ;
_weightVB[BB][6]  =  0.30864197530864    ;
_weightVB[BB][7]  =  0.49382716049383    ;
_weightVB[BB][8]  =  0.30864197530864    ;
                                         
                                      
		      break;
		    }

		    case(2) : {

// ========= VV ========
_weightVB[VV][0] =   0.30864197530864    ;
_weightVB[VV][1] =   0.49382716049383    ;
_weightVB[VV][2] =   0.30864197530864    ;
_weightVB[VV][3] =   0.49382716049383    ;
_weightVB[VV][4] =   0.79012345679012    ;
_weightVB[VV][5] =   0.49382716049383    ;
_weightVB[VV][6]  =  0.30864197530864    ;
_weightVB[VV][7]  =  0.49382716049383    ;
_weightVB[VV][8]  =  0.30864197530864    ;

// ========= BB ========
_weightVB[BB][0] =   0.555555555555555555555555555556    ;
_weightVB[BB][1] =   0.888888888888888888888888888889    ;
_weightVB[BB][2] =   0.555555555555555555555555555556    ;

		      break;
		    }

		    default: {std::cout << "FE: Dimension not implemented" << std::endl; abort();}
		    
		  } //end dim
	   
	   break;
         } //end case quadrilateral
      
	 case(TRIANG) :  {
               switch(_geomel->_dim) {
		    
		    case(3) : {
		      
// ========= VV ========  //these are not even the same number
_weightVB[VV][0] =    0.030283678097089;
_weightVB[VV][1] =    0.006026785714286;  
_weightVB[VV][2] =    0.006026785714286;  
_weightVB[VV][3] =    0.006026785714286;
_weightVB[VV][4] =    0.006026785714286;
_weightVB[VV][5] =    0.011645249086029;
_weightVB[VV][6]  =   0.011645249086029;
_weightVB[VV][7]  =   0.011645249086029;
_weightVB[VV][8]  =   0.011645249086029;
_weightVB[VV][9]  =   0.010949141561386;
_weightVB[VV][10] =   0.010949141561386;
_weightVB[VV][11] =   0.010949141561386;
_weightVB[VV][12] =   0.010949141561386;
_weightVB[VV][13] =   0.010949141561386;     
_weightVB[VV][14] =   0.010949141561386;
		      
		      
// ========= BB ========  //these are multiplied by 4 and not with the same order
_weightVB[BB][0] =   0.1125;
_weightVB[BB][1] =   0.062969590272414;
_weightVB[BB][2] =   0.062969590272414;
_weightVB[BB][3] =   0.062969590272414;
_weightVB[BB][4] =   0.066197076394253;
_weightVB[BB][5] =   0.066197076394253;
_weightVB[BB][6]  =  0.066197076394253;
		      
       
		      break;
		    }

		    case(2) : {
		      
// ========= VV ========   //these are multiplied by 4 and not with the same order
_weightVB[VV][0] =   0.1125;  
_weightVB[VV][1] =   0.062969590272414;  
_weightVB[VV][2] =   0.062969590272414;  
_weightVB[VV][3] =   0.062969590272414;  
_weightVB[VV][4] =   0.066197076394253;  
_weightVB[VV][5] =   0.066197076394253;  
_weightVB[VV][6]  =  0.066197076394253;  

// ========= BB ========
_weightVB[BB][0] =   0.555555555555555555555555555556    ;
_weightVB[BB][1] =   0.888888888888888888888888888889    ;
_weightVB[BB][2] =   0.555555555555555555555555555556    ;

		      break;
		    }

		    default: {std::cout << "FE: Dimension not implemented" << std::endl; abort();}
		    
		  } //end dim 

	   break;
          } //end case triangular
      
       } //end switch geomel   
  
  
       
 }
 
 
 


} //end namespace femus


 