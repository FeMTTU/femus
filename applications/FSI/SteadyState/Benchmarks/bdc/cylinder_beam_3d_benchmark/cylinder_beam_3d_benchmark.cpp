#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>


static const double inflow_elong = 0.3;
static const double um = 0.2;
static const double L = 2.5 + inflow_elong; //elongated
static const double H = 0.41;
static const double factor_avg = /*1.5*/36. /*9/4 * 16*/;



extern "C" double InitalValueU(const std::vector < double >& x) {
    
  double xc = 0.2 + inflow_elong;
  double yc = 0.2;
  double r = 0.05;
  double r2 = r * r;
  double xMxc2 = (x[0] - xc) * (x[0] - xc);
  double OMxc2 = (0. - xc) * (0. - xc);
  double yMyc2 = (x[1] - yc) * (x[1] - yc);

  return (xMxc2 + yMyc2 - r2)/(OMxc2 + yMyc2 - r2) * ( factor_avg * um * 1./(H * H * H * H) * x[1] * (H - x[1]) ) * x[2] * (H - x[2]) * exp( - L * x[0] );

}

//---------------------------------------------------------------------------------------------------------------------

extern "C" bool BdcFunction(const std::vector < double >& x,const char name[], double &value, const int facename, const double time) {
    
  bool test = 1; //dirichlet
  value = 0.;
  if(!strcmp(name,"U")) {
      
    if(1 == facename){   //inflow
      test = 1;
      value = factor_avg * um * 1. /(H * H * H * H) * x[1] * (H - x[1]) * x[2] * (H - x[2]);

    }
    else if(2 == facename ){  //outflow
      test = 0;  // u \cdot n free
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
      test = 0;   // u \times n free
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
  else if(!strcmp(name,"W")) {
      
    if( 1 == facename){            //inflow
      test = 1;
      value = 0.;
    }
    else if( 2 == facename ){      //outflow
      test = 0;  // u \times n free
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
  else if(!strcmp(name,"P")) {
      
    if(1 == facename){
      test = 0;
      value = 0.;
    }
    else if( 2 == facename ){      //outflow
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
  else if(!strcmp(name,"DX")) {
      
    if( 1 == facename ){         //inflow
      test = 1;
      value = 0.;
    }
    else if( 2 == facename ){   //outflow
      test = 1;
      value = 0.;
    }
    else if( 3 == facename ){   // no-slip fluid wall
      test = 0; ///@todo do I have to leave the displacement of the fluid here free?
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
  else if(!strcmp(name,"DY")) {
      
    if( 1 == facename ){         //inflow
      test = 0; ///@todo do I have to leave the displacement of the fluid here free?
      value = 0.;
    }
    else if( 2 == facename ){   //outflow
      test = 0; ///@todo do I have to leave the displacement of the fluid here free?
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
  else if(!strcmp(name,"DZ")) {
      
    if( 1 == facename ){         //inflow
      test = 0; ///@todo do I have to leave the displacement of the fluid here free?
      value = 0.;
    }
    else if( 2 == facename ){   //outflow
      test = 0; ///@todo do I have to leave the displacement of the fluid here free?
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
