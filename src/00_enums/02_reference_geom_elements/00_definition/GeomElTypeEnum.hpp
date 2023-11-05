#ifndef __femus_enums_GeomElTypeEnum_hpp__
#define __femus_enums_GeomElTypeEnum_hpp__

#include <string>
#include <vector>

namespace femus {
    
  static const std::vector< std::string >  geom_elems = {"hex", "tet", "wedge", "quad", "tri", "line"};

}

enum  GeomElType { HEX=0, 
                   TET,   //1
		   WEDGE, //2
		   QUAD,  //3
		   TRI,   //4
		   LINE   //5
  
                 };

                 
                
#define N_GEOM_ELS 6

#define MAXIMUM_NUMBER_OF_FACES_PER_GEOM_EL   6
#define MAXIMUM_NUMBER_OF_NODES_PER_FACE      9

#define MAXIMUM_NUMBER_OF_NON_BIQUADRATIC_NODES_TO_USE_IN_WEIGHTING   18

#define MAX_EL_N_NODES 27

#define MAX_EL_N_FACES 6
		 
#define LEVEL_AT_WHICH_YOU_PICK_THE_DIM  0

#define ZERO_ELEM   0
		 
#define ZERO_FACE   0

#endif
