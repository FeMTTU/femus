#include "GeomEl.hpp"
#include "GeomElTypeEnum.hpp"
#include "Elem.hpp"

#include <iostream>
#include <cstdlib>
#include <cstring>


namespace femus {

// NRE[6]= {8,8,8,4,4,2};

  GeomEl::GeomEl(const std::string geomel_id_in, const uint mesh_order) :
     _geomel_id(geomel_id_in)    {

     if (mesh_order != QQ && mesh_order != LL ) { std::cout << "Wrong mesh order" << std::endl; abort(); }  
     
     if (!strcmp(_geomel_id.c_str(),"hex")) {  
           n_se = NRE[HEX];
      switch(mesh_order) {
      case(QQ):
          _elnds = 27; break;
      case(LL):
          _elnds = 8;  break;
        } 
      } //end hex

      else if (!strcmp(_geomel_id.c_str(),"wedge")) {
           std::cout << "Not supported yet" << std::endl; abort();
      } //end wedge
      
      else if (!strcmp(_geomel_id.c_str(),"tet")) {
           n_se = NRE[TET];
      switch(mesh_order) {
      case(QQ):
          _elnds = 10; break;
      case(LL):
          _elnds = 4;  break;
        } 
      } //end tet
      
      else if (!strcmp(_geomel_id.c_str(),"quad")) {
           n_se = NRE[QUAD];
      switch(mesh_order) {
      case(QQ):
          _elnds = 9; break;
      case(LL):
          _elnds = 4; break;
        } 
      }  //end quad
      
      else if (!strcmp(_geomel_id.c_str(),"tri")) {
           n_se = NRE[TRI];
      switch(mesh_order) {
      case(QQ):
          _elnds = 6; break;
      case(LL):
          _elnds = 3; break;
        } 
      }  //end tri
      
      else if (!strcmp(_geomel_id.c_str(),"line")) { 
           n_se = NRE[LINE];
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
