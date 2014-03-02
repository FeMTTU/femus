#include "QRule.hpp"

#include <cstdlib>

 
#include "GeomElTypeEnum.hpp"
#include "GeomEl.hpp"
 

  QRule::~QRule()  {
    
 for (uint vb=0; vb < VB; vb++)  {
    delete []        _weightVB[vb];
  }
  
 }
 
 
 QRule::QRule(GeomEl* geomel_in)  : _geomel(geomel_in)  {
   
   
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
		      _NoGaussVB[VV] = 14;
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
_weightVB[VV][0] =   0.17146776410151     ; 
_weightVB[VV][1] =   0.27434842257476     ;   
_weightVB[VV][2] =   0.17146776410151     ;   
_weightVB[VV][3] =     0.27434842257476    ;
_weightVB[VV][4] =     0.43895747613937    ;
_weightVB[VV][5] =     0.27434842257476    ;
_weightVB[VV][6]  =    0.17146776410151    ;
_weightVB[VV][7]  =    0.27434842257476    ;
_weightVB[VV][8]  =    0.17146776410151    ;
_weightVB[VV][9]  =    0.27434842257476    ;
_weightVB[VV][10] =    0.43895747613937    ;
_weightVB[VV][11] =     0.27434842257476   ;
_weightVB[VV][12] =     0.43895747613937   ;
_weightVB[VV][13] =     0.7023319618546    ;
_weightVB[VV][14] =     0.43895747613937   ;
_weightVB[VV][15] =     0.27434842257476   ;
_weightVB[VV][16] =     0.43895747613937   ;
_weightVB[VV][17] =     0.27434842257476   ;
_weightVB[VV][18] =     0.17146776410151   ;
_weightVB[VV][19] =     0.27434842257476   ;
_weightVB[VV][20] =     0.17146776410151   ;
_weightVB[VV][21] =     0.27434842257476   ;
_weightVB[VV][22] =     0.43895747613937   ;
_weightVB[VV][23] =     0.27434842257476   ;
_weightVB[VV][24] =     0.17146776410151   ;
_weightVB[VV][25] =     0.27434842257476   ;
_weightVB[VV][26] =     0.17146776410151   ;
		      
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
_weightVB[BB][1] =   0.555555555555555555555555555556    ;
_weightVB[BB][2] =   0.888888888888888888888888888889    ;

		      break;
		    }

		    default: {std::cout << "FE: Dimension not implemented" << std::endl; abort();}
		    
		  } //end dim
	   
	   break;
         } //end case quadrilateral
      
	 case(TRIANG) :  {
               switch(_geomel->_dim) {
		    
		    case(3) : {
		      
// ========= VV ========
_weightVB[VV][0] =   0.0375626419060053  ; 
_weightVB[VV][1] =   0.0375626419060053  ;   
_weightVB[VV][2] =   0.0375626419060053  ;   
_weightVB[VV][3] =   0.0375626419060053  ;
_weightVB[VV][4] =   0.0244976810387873  ;
_weightVB[VV][5] =   0.0244976810387873  ;
_weightVB[VV][6]  =  0.0244976810387873  ;
_weightVB[VV][7]  =  0.0244976810387873  ;
_weightVB[VV][8]  =  0.0141820069256938  ;
_weightVB[VV][9]  =  0.0141820069256938  ;
_weightVB[VV][10] =  0.0141820069256938  ;
_weightVB[VV][11] =  0.0141820069256938  ;
_weightVB[VV][12] =  0.0141820069256938  ;
_weightVB[VV][13] =  0.0141820069256938  ;      
		      
		      
// ========= BB ========
_weightVB[BB][0] =   0.45        ;
_weightVB[BB][1] =   0.264788305577012        ;
_weightVB[BB][2] =   0.264788305577012        ;
_weightVB[BB][3] =   0.264788305577012        ;
_weightVB[BB][4] =   0.251878361089654        ;
_weightVB[BB][5] =   0.251878361089654        ;
_weightVB[BB][6]  =  0.251878361089654        ;
		      
      
		      break;
		    }

		    case(2) : {
		      
// ========= VV ========
_weightVB[VV][0] =   0.45        ;
_weightVB[VV][1] =   0.264788305577012        ;
_weightVB[VV][2] =   0.264788305577012        ;
_weightVB[VV][3] =   0.264788305577012        ;
_weightVB[VV][4] =   0.251878361089654        ;
_weightVB[VV][5] =   0.251878361089654        ;
_weightVB[VV][6]  =  0.251878361089654        ;

// ========= BB ========
_weightVB[BB][0] =   0.555555555555555555555555555556    ;
_weightVB[BB][1] =   0.555555555555555555555555555556    ;
_weightVB[BB][2] =   0.888888888888888888888888888889    ;

		      break;
		    }

		    default: {std::cout << "FE: Dimension not implemented" << std::endl; abort();}
		    
		  } //end dim 

	   break;
          } //end case triangular
      
       } //end switch geomel 
  
      
       
       
 }
 
 
 
 