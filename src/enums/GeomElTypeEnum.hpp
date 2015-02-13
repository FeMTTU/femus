#ifndef _geomeltype_enum
#define _geomeltype_enum


enum  GeomElType { HEX=0, 
                   TET,   //1
		   WEDGE, //2
		   QUAD,  //3
		   TRI,   //4
		   LINE   //5
  
                 };

#define N_GEOM_ELS 6
		 
#define LEV_PICK  0


#endif