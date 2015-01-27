#include "GeomEl.hpp"
#include "GeomElTypeEnum.hpp"
#include "Elem.hpp"

#include <iostream>
#include <cstdlib>
#include <cstring>


namespace femus {  
  //NVE[][2] for biquadratic
  //NVE[][0] for linear

// // //vertices,edges,faces,interior,element,element+derivatives
// //   const unsigned NVE[6][5]= {
// //     {8,20,27,1,4},  //hex
// //     {4,10,10,1,4},   //tet
// //     {6,15,18,1,4},   //wedge
// //     {4, 8, 9,1,3},   //quad
// //     {3, 6, 6,1,3},   //tri
// //     {2, 3, 3,1,2}    //line
// //   };

  GeomEl::GeomEl(const std::string geomel_id_in, const uint mesh_order) {

//      if (mesh_order != QQ && mesh_order != LL ) { std::cout << "Wrong mesh order" << std::endl; abort(); }  
//      
//      if (!strcmp(geomel_id_in.c_str(),"hex")) {  
//       switch(mesh_order) {
//       case(QQ):
//       case(LL):
//         } 
//       } //end hex
// 
//       else if (!strcmp(geomel_id_in.c_str(),"wedge")) {
//            std::cout << "Not supported yet" << std::endl; abort();
//       } //end wedge
//       
//       else if (!strcmp(geomel_id_in.c_str(),"tet")) {
//       switch(mesh_order) {
//       case(QQ):
//       case(LL):
//         } 
//       } //end tet
//       
//       else if (!strcmp(geomel_id_in.c_str(),"quad")) {
//       switch(mesh_order) {
//       case(QQ):
//       case(LL):
//         } 
//       }  //end quad
//       
//       else if (!strcmp(geomel_id_in.c_str(),"tri")) {
//       switch(mesh_order) {
//       case(QQ):
//       case(LL):
//         } 
//       }  //end tri
//       
//       else if (!strcmp(geomel_id_in.c_str(),"line")) { 
//       switch(mesh_order) {
//       case(QQ):
//       case(LL):
//         } 
//       }  //end line
//       
//       else {std::cout << "Geometric element not recognized\n"; abort(); }


  } 

  
  GeomEl::~GeomEl() {   }


} //end namespace femus
