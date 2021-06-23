#ifndef __femus_enums_GeomElTypeEnum_hpp__
#define __femus_enums_GeomElTypeEnum_hpp__

#include <string>
namespace femus {
    
  static const std::vector< std::string >  geom_elems = {"hex", "tet", "wedge", "quad", "tri", "line"};
  static const std::vector< std::string >  fe_fams    = {"linear", "quadratic", "biquadratic", "constant", "disc_linear"};

}

enum  GeomElType { HEX=0, 
                   TET,   //1
		   WEDGE, //2
		   QUAD,  //3
		   TRI,   //4
		   LINE   //5
  
                 };

#define N_GEOM_ELS 6

#define MAX_EL_N_NODES 27

#define MAX_EL_N_FACES 6
		 
#define LEV_PICK  0

#define ZERO_ELEM   0
		 
#define ZERO_FACE   0

#endif
