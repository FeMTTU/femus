#include "GeomEl.hpp"
#include "GeomElTypeEnum.hpp"
#include "Elem.hpp"

#include <iostream>
#include <cstdlib>
#include <cstring>


namespace femus {


  GeomEl::GeomEl(const std::string geomel_id_in, const uint mesh_order) {

     if (mesh_order != QQ && mesh_order != LL ) { std::cout << "Wrong mesh order" << std::endl; abort(); }  
     
     if (!strcmp(geomel_id_in.c_str(),"hex")) {  
      switch(mesh_order) {
      case(QQ):
          _elnds = 27; break;
      case(LL):
          _elnds = 8;  break;
        } 
      } //end hex

      else if (!strcmp(geomel_id_in.c_str(),"wedge")) {
           std::cout << "Not supported yet" << std::endl; abort();
      } //end wedge
      
      else if (!strcmp(geomel_id_in.c_str(),"tet")) {
      switch(mesh_order) {
      case(QQ):
          _elnds = 10; break;
      case(LL):
          _elnds = 4;  break;
        } 
      } //end tet
      
      else if (!strcmp(geomel_id_in.c_str(),"quad")) {
      switch(mesh_order) {
      case(QQ):
          _elnds = 9; break;
      case(LL):
          _elnds = 4; break;
        } 
      }  //end quad
      
      else if (!strcmp(geomel_id_in.c_str(),"tri")) {
      switch(mesh_order) {
      case(QQ):
          _elnds = 6; break;
      case(LL):
          _elnds = 3; break;
        } 
      }  //end tri
      
      else if (!strcmp(geomel_id_in.c_str(),"line")) { 
      switch(mesh_order) {
      case(QQ):
          _elnds = 3; break;
      case(LL):
          _elnds = 2; break;
        } 
      }  //end line
      
      else {std::cout << "Geometric element not recognized\n"; abort(); }


  } 

  
  GeomEl::~GeomEl() {   }


} //end namespace femus
