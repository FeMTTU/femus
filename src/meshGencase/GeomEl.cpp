#include "GeomEl.hpp"
#include "GeomElTypeEnum.hpp"

#include <iostream>
#include <cstdlib>
#include <cstring>


namespace femus {



//This class was badly constructed 
//because there were things only quadratic and things both quadratic and linear

  GeomEl::GeomEl(const std::string geomel_id_in, const uint mesh_order) :
     _geomel_id(geomel_id_in)    {

     if (mesh_order != QQ && mesh_order != LL ) { std::cout << "Wrong mesh order" << std::endl; abort(); }  
     
     if (!strcmp(_geomel_id.c_str(),"hex")) {  
           _dim = 3;
           n_se = 8;
      switch(mesh_order) {
      case(QQ):
          _elnds = 27; name = "Hexahedron_27";  break;
      case(LL):
          _elnds = 8;  name = "Hexahedron";     break;
        } 
      } //end hex

      else if (!strcmp(_geomel_id.c_str(),"wedge")) {
           std::cout << "Not supported yet" << std::endl; abort();
      } //end wedge
      
      else if (!strcmp(_geomel_id.c_str(),"tet")) {
           _dim = 3;
           n_se = 8;
      switch(mesh_order) {
      case(QQ):
          _elnds = 10; name = "Tetrahedron_10";  break;
      case(LL):
          _elnds = 4;  name = "Tetrahedron";     break;
        } 
      } //end tet
      
      else if (!strcmp(_geomel_id.c_str(),"quad")) {
           _dim = 2;
           n_se = 4;
      switch(mesh_order) {
      case(QQ):
          _elnds = 9; name = "Quadrilateral_9";  break;
      case(LL):
          _elnds = 4; name = "Quadrilateral";     break;
        } 
      }  //end quad
      
      else if (!strcmp(_geomel_id.c_str(),"tri")) {
           _dim = 2;
           n_se = 4;
      switch(mesh_order) {
      case(QQ):
          _elnds = 6; name = "Triangle_6";  break;
      case(LL):
          _elnds = 3; name = "Triangle";     break;
        } 
      }  //end tri
      
      else if (!strcmp(_geomel_id.c_str(),"line")) { 
           _dim = 1;
           n_se = 2;
      switch(mesh_order) {
      case(QQ):
          _elnds = 3; name = "Edge_3";  break;
      case(LL):
          _elnds = 2; name = "Edge";     break;
        } 
      }  //end line
      
      else {std::cout << "Geometric element not recognized\n"; abort(); }


  } 
   //====end constructor
   
  GeomEl::~GeomEl() {   }
  //====end destructor


} //end namespace femus


// #endif
// #endif


   