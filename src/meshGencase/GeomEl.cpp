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
           _dim = 3;
           n_se = NRE[HEX];
      switch(mesh_order) {
      case(QQ):
          _elnds = 27; _xdmf_name = "Hexahedron_27";  break;
      case(LL):
          _elnds = 8;  _xdmf_name = "Hexahedron";     break;
        } 
      } //end hex

      else if (!strcmp(_geomel_id.c_str(),"wedge")) {
           std::cout << "Not supported yet" << std::endl; abort();
      } //end wedge
      
      else if (!strcmp(_geomel_id.c_str(),"tet")) {
           _dim = 3;
           n_se = NRE[TET];
      switch(mesh_order) {
      case(QQ):
          _elnds = 10; _xdmf_name = "Tetrahedron_10";  break;
      case(LL):
          _elnds = 4;  _xdmf_name = "Tetrahedron";     break;
        } 
      } //end tet
      
      else if (!strcmp(_geomel_id.c_str(),"quad")) {
           _dim = 2;
           n_se = NRE[QUAD];
      switch(mesh_order) {
      case(QQ):
          _elnds = 9; _xdmf_name = "Quadrilateral_9";  break;
      case(LL):
          _elnds = 4; _xdmf_name = "Quadrilateral";     break;
        } 
      }  //end quad
      
      else if (!strcmp(_geomel_id.c_str(),"tri")) {
           _dim = 2;
           n_se = NRE[TRI];
      switch(mesh_order) {
      case(QQ):
          _elnds = 6; _xdmf_name = "Triangle_6";  break;
      case(LL):
          _elnds = 3; _xdmf_name = "Triangle";     break;
        } 
      }  //end tri
      
      else if (!strcmp(_geomel_id.c_str(),"line")) { 
           _dim = 1;
           n_se = NRE[LINE];
      switch(mesh_order) {
      case(QQ):
          _elnds = 3; _xdmf_name = "Edge_3";  break;
      case(LL):
          _elnds = 2; _xdmf_name = "Edge";     break;
        } 
      }  //end line
      
      else {std::cout << "Geometric element not recognized\n"; abort(); }


  } 

  
  GeomEl::~GeomEl() {   }


} //end namespace femus
