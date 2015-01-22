//C++
#include <cmath>

//library headers
#include "Box.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "MultiLevelMeshTwo.hpp" 
#include "NormTangEnum.hpp"
#include "TimeLoop.hpp"

//application headers
#include "TempQuantities.hpp"
#include "EqnT.hpp"


// In these files, just like in the equation file, you need to have an EXPLICIT KNOWLEDGE 
// of the ORDER of the UNKNOWNS

//If for a given equation you decide to change the NUMBER of UNKNOWNS then it is a little bit of a problem
//in the sense that you cannot do it totally automatically, you need to revisit the EQUATION and the BC/IC.


void EqnT::ic_read(const double xp[], double u_value[], const double el_xm[]) const {

  Box* box = static_cast<Box*>(_mesh.GetDomain());
  
  double* x_rotshift = new double[_mesh.get_dim()];
  _mesh._domain->TransformPointToRef(xp,x_rotshift); 

    u_value[0] = 3.; 
    u_value[1] = 4.;  
    u_value[2] = 5.;

  delete[] x_rotshift;
  
  return;
  
}



void EqnT::bc_read(const double xp[],const double /*normal */[],int bc_flag[]) const {

  const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");
  

   Box* box = static_cast<Box*>(_mesh.GetDomain());

  double* lb = new double[_mesh.get_dim()];
  double* le = new double[_mesh.get_dim()];
  lb[0] = box->_lb[0];  //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
   if ( _mesh.get_dim() == 3 ) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
   }
  
  double* x_rotshift = new double[_mesh.get_dim()];
  _mesh._domain->TransformPointToRef(xp,x_rotshift);
  
   if ( _mesh.get_dim() == 2 ) {
  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {  //left of the RefBox

  bc_flag[0]=0;
  bc_flag[1]=0;  
  bc_flag[2]=0;

  }

 if ( (le[0]-lb[0])  - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll)  { //right of the RefBox

   bc_flag[0]=0;
   bc_flag[1]=0; 
   bc_flag[2]=0;
   
   }

  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox

//       bc_flag[0]=0; 
//       bc_flag[1]=0; 
//       bc_flag[2]=0; 

  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox

//       bc_flag[0]=0;
//       bc_flag[1]=0; 
//       bc_flag[2]=0;
      
  }

   } //dim 2

  else if ( _mesh.get_dim() == 3 ) {
    
  
  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox

    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0; 
    
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox

    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
  
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
//      bc_flag[0]=0;    
//      bc_flag[1]=0; 
//      bc_flag[2]=0; 

  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
//      bc_flag[0]=0;
//      bc_flag[1]=0;
//      bc_flag[2]=0;
     
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
//      bc_flag[0]=0;
//      bc_flag[1]=0;
//      bc_flag[2]=0;

  }
 
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
//       bc_flag[0]=0;
//       bc_flag[1]=0;
//       bc_flag[2]=0;

  }
  
  }  // dim 3
  

  delete[] lb;
  delete[] le;
  delete[] x_rotshift;
  
  return;
}
