#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>

// double InitalValueU(const std::vector < double >& x);
//
// bool SetBoundaryConditionTurek_2D_FSI_and_solid(const std::vector < double >& x,const char name[],
// 						double &value, const int FaceName, const double = 0.);

const double um = 0.2;
const double L = 2.5;
const double H = 0.41;

extern "C" double InitalValueU(const std::vector < double >& x) {
  double xc = 0.2;
  double yc = 0.2;
  double r = 0.05;
  double r2 = r * r;
  double xMxc2 = (x[0] - xc) * (x[0] - xc);
  double OMxc2 = (0. - xc) * (0. - xc);
  double yMyc2 = (x[1] - yc) * (x[1] - yc);

  return (xMxc2 + yMyc2 - r2)/(OMxc2 + yMyc2 - r2) * ( 1.5 * um * 4.0/(H * H) * x[1] * (H - x[1]) ) * exp( - L * x[0] );

}

//---------------------------------------------------------------------------------------------------------------------

extern "C" bool BdcFunction(const std::vector < double >& x,const char name[], double &value, const int facename, const double time) {

  bool test = 1; //dirichlet
  value = 0.;
  
  
  if(!strcmp(name,"U")) {
      
    if(1 == facename){   //inflow
      test = 1;
      value = 1.5 * um * 4.0/(H * H) * x[1] * (H - x[1]);

    }
    else if(2 == facename ){  //outflow
      test = 0;
      //    test=1;
      value=0.;
    }
    else if(3 == facename ){  // no-slip fluid wall
      test = 1;
      value = 0.;
    }
    else if(4 == facename ){  // no-slip solid wall
      test = 1;
      value = 0.;
    }
    else if( 6 == facename ){   // beam case zero stress
      test = 0;
      value = 0.;
    }
    
  }
  else if(!strcmp(name,"V")) {
      
    if( 1 == facename){            //inflow
      test = 1;
      value = 0.;
    }
    else if( 2 == facename ){      //outflow
      test = 0;
      //    test=1;
      value = 0.;
    }
    else if( 3 == facename ){      // no-slip fluid wall
      test = 1;
      value = 0;
    }
    else if( 4 == facename ){      // no-slip solid wall
      test = 1;
      value = 0.;
    }
    else if( 6 == facename ){   // beam case zero stress
      test = 0;
      value = 0.;
    }
    
  }
  else if(!strcmp(name,"P")){
      
    if(1 == facename){
      test = 0;
      value = 0.;
    }
    else if( 2 == facename ){
      test = 0;
      value = 0.;
    }
    else if( 3 == facename ){
      test = 0;
      value = 0.;
    }
    else if( 4 == facename ){
      test = 0;
      value = 0.;
    }
    else if( 6 == facename ){   // beam case zero stress
      test = 0;
      value = 0.;
    }
    
  }
  else if(!strcmp(name,"DX")){
      
    if( 1 == facename ){         //inflow
      test = 1;
      value = 0.;
    }
    else if( 2 == facename ){   //outflow
      test = 1;
      value = 0.;
    }
    else if( 3 == facename ){   // no-slip fluid wall
      test = 0; //0
      value = 0.;
    }
    else if( 4 == facename ){   // no-slip solid wall
      test = 1;
      value = 0.;
    }
    else if( 6 == facename ){   // beam case zero stress
      test = 0;
      value = 0.;
    }
    
  }
  else if(!strcmp(name,"DY")){
      
    if( 1 == facename ){         //inflow
      test = 0; // 0
      value = 0.;
    }
    else if( 2 == facename ){   //outflow
      test = 0; // 0
      value = 0.;
    }
    else if( 3 == facename ){   // no-slip fluid wall
      test = 1;
      value = 0.;
    }
    else if( 4 == facename ){   // no-slip solid wall
      test = 1;
      value = 0.;
    }
    else if( 6 == facename ){   // beam case zero stress
      test = 0;
      value = 0.;
    }
    
  }
  return test;
}
