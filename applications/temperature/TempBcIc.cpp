//C++
#include <cmath>

//library headers
#include "Box.hpp"
#include "EquationsMap.hpp"
#include "Physics.hpp"
#include "MeshTwo.hpp" 
#include "NormTangEnum.hpp"
#include "TimeLoop.hpp"

//application headers
#include "Temp_conf.hpp"
#include "TempQuantities.hpp"
#include "TempPhysics.hpp"
#include "EqnT.hpp"


void EqnT::ic_read(double xp[],double u_value[], double el_xm[]) {

//     const double Tref = _phys._physrtmap.get("Tref");
  const double bdry_toll = _mesh._mesh_rtmap.get("bdry_toll");


  Box* box= static_cast<Box*>(_mesh.GetDomain());
  
  double* x_rotshift = new double[_mesh.get_dim()];
  _mesh._domain->TransformPointToRef(xp,x_rotshift); 

#if (DIMENSION==2)

//   const double magnitude = (x_rotshift[0] - box->_lb[0])*(box->_le[0]-x_rotshift[0])*(x_rotshift[1] - box->_lb[1])/(Tref); 

/*T'*/    u_value[0] = 0.; 
/*T_0*/   u_value[1] = 1.;  //0.*magnitude
  if ((box->_le[1]-box->_lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (box->_le[1]-box->_lb[1]) -(x_rotshift[1]) < bdry_toll)  { u_value[1] = 0.; }   

/*T_adj*/ u_value[2] = 0.;
#if FOURTH_ROW==1
/*p2*/    u_value[3] =  72.*(xp[0]);
#endif


#elif (DIMENSION==3)  
    
    u_value[0] = 0.;
    u_value[1] = 0.;
    u_value[2] = 0.;
#if FOURTH_ROW==1
    u_value[3] =  72.*(xp[0]);
#endif
    
    
#endif    
    
  delete[] x_rotshift;
  
  return;
}



void EqnT::bc_read(double xp[],double /*normal */[],int bc_flag[]) {
// T' and its adjoint must be Dirichlet homogeneous everywhere on the boundary, by definition.


  const double bdry_toll = _mesh._mesh_rtmap.get("bdry_toll");
  

Box* box= static_cast<Box*>(_mesh.GetDomain());

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
  
  
#if (DIMENSION==2)
  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {  //left of the RefBox

  bc_flag[0]=0;  //always fixed //T'
  
    if  (_eqnmap._timeloop._curr_t_idx < 1)  
   {   bc_flag[1]=0; }
    else if ( (x_rotshift[1]) < 0.4*(le[1]-lb[1])  ||  (x_rotshift[1]) > 0.6*(le[1]-lb[1]) )  {  bc_flag[1]=0;  } 
    //T_0

  bc_flag[2]=0;  //always fixed//T_adj
#if FOURTH_ROW==1
   bc_flag[3]=0;
#endif  
   ///////////// ////////////////////
//    if ( (x_rotshift[1]) < 0.25*(le[1]-lb[1])  ||  (x_rotshift[1]) > 0.75*(le[1]-lb[1]) )  {  bc_flag[1]=0;  } 
   ///////////////////////////////////
 
  }

 if ( (le[0]-lb[0])  - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll)  { //right of the RefBox

//    if ( (x_rotshift[1]) < 0.4*(le[1]-lb[1])  ||  (x_rotshift[1]) > 0.6*(le[1]-lb[1]) )  { 
   bc_flag[0]=0;
   bc_flag[1]=0; 
   bc_flag[2]=0;
//     } 
   
   }

  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox

      bc_flag[0]=0; //always fixed

      //===== START FROM A SIMPLE STATE SOLUTION =========
if  (_eqnmap._timeloop._curr_t_idx < 1)   bc_flag[1]=0;     //      bc_flag[1]=0; //=====CONTROL//////////////
      //===== START FROM A SIMPLE STATE SOLUTION =========
else if  ( (x_rotshift[0]) < 0.25*(le[0] - lb[0]) || ( x_rotshift[0]) > 0.75*(le[0] - lb[0]) )  {  bc_flag[1]=0; }
      
      bc_flag[2]=0; //always fixed

  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox

      bc_flag[0]=0; //always fixed

      bc_flag[1]=0; 

//  if  (  _eqnmap._timeloop._curr_t_idx < 1 )
// {  bc_flag[1]=0; }
//  else {
//    if ( (x_rotshift[0]) <= 0.70*(le[0]-lb[0]) )  {  bc_flag[1]=0; }
//     else {  if (bc_flag[1]==0) bc_flag[1]=1;  } //se qualcuno ha messo uno zero li' mettici 1, voglio far vincere il controllo!
//  }
  
      bc_flag[2]=0; //always fixed
      
  }


  
#elif (DIMENSION==3)

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox

    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0; 
    
#if FOURTH_ROW==1
        bc_flag[3]=0;
#endif    
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox

    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;

//         bc_flag[3]=0;
    
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
     bc_flag[0]=0;    
     bc_flag[1]=0; 
     bc_flag[2]=0; 

//              bc_flag[3]=0;
	     
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
     bc_flag[0]=0;
     bc_flag[1]=0;
     bc_flag[2]=0;
     
//              bc_flag[3]=0;
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
     bc_flag[0]=0;
     bc_flag[1]=0;
     bc_flag[2]=0;
//              bc_flag[3]=0;
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
      bc_flag[0]=0;
      bc_flag[1]=0;
      bc_flag[2]=0;
//             bc_flag[3]=0;
  }
  
  
#endif
  
  delete[] lb;
  delete[] le;
  delete[] x_rotshift;
  
  return;
}
