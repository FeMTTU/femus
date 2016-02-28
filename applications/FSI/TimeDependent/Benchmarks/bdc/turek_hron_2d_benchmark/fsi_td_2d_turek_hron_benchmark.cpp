#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>

const double um = 2.0;
const double L = 2.5;
const double H = 0.41;

//---------------------------------------------------------------------------------------------------------------------

extern "C" bool BdcFunction(const std::vector < double >& x,const char name[], double &value, const int facename, const double time) {
  bool test=1; //dirichlet
  value=0.;
  if( !strcmp(name,"U") ) {
    if(1==facename){   //inflow
      test = 1;
      if(time < 2.0) {
        value = 1.5 * um * 4.0 / 0.1681 * x[1] * ( H - x[1]) * 0.5 * ( 1. - cos( 0.5 * 3.141592653589793 * time) );
      }
      else {
        value = 1.5 * um * 4.0 / 0.1681 * x[1] * ( H - x[1] );
      }
    }
    else if(2==facename ){  //outflow
      test=0;
      //    test=1;
      value=0.;
    }
    else if(3==facename ){  // no-slip fluid wall
      test=1;
      value=0.;
    }
    else if(4==facename ){  // no-slip solid wall
      test=1;
      value=0.;
    }
    else if(6==facename ){   // beam case zero stress
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"V")){
    if(1==facename){            //inflow
      test=1;
      value=0.;
    }
    else if(2==facename ){      //outflow
      test=0;
      //    test=1;
      value=0.;
    }
    else if(3==facename ){      // no-slip fluid wall
      test=1;
      value=0;
    }
    else if(4==facename ){      // no-slip solid wall
      test=1;
      value=0.;
    }
    else if(6==facename ){   // beam case zero stress
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"P")){
    if(1==facename){
      test=0;
      value=0.;
    }
    else if(2==facename ){
      test=0;
      value=0.;
    }
    else if(3==facename ){
      test=0;
      value=0.;
    }
    else if(4==facename ){
      test=0;
      value=0.;
    }
    else if(6==facename ){   // beam case zero stress
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"DX")){
    if(1==facename){         //inflow
      test=1;
      value=0.;
    }
    else if(2==facename ){   //outflow
      test=1;
      value=0.;
    }
    else if(3==facename ){   // no-slip fluid wall
      test=0;
      value=0.;
    }
    else if(4==facename ){   // no-slip solid wall
      test=1;
      value=0.;
    }
    else if(6==facename ){   // beam case zero stress
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"DY")){
    if(1==facename){         //inflow
      test=0; // 0
      value=0.;
    }
    else if(2==facename ){   //outflow
      test=0; // 0
      value=0.;
    }
    else if(3==facename ){   // no-slip fluid wall
      test=1;
      value=0.;
    }
    else if(4==facename ){   // no-slip solid wall
      test=1;
      value=0.;
    }
    else if(6==facename ){   // beam case zero stress
      test=0;
      value=0.;
    }
  }
  return test;
}