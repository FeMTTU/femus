#include <cmath>

//library headers
#include "Box.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "MultiLevelMeshTwo.hpp" 
#include "NormTangEnum.hpp"


//application headers
#include "Opt_conf.hpp"
#include "OptQuantities.hpp"
#include "EqnNS.hpp"
#include "EqnNSAD.hpp"
#include "EqnMHD.hpp"
#include "EqnMHDAD.hpp"
#include "EqnMHDCONT.hpp"

namespace femus {


/// This function generates the initial conditions for the NS system:
//xp[] is the NON-DIMENSIONAL node coordinate
//when this function is called,
//the domain has been NONDIMENSIONALIZED
//but NOT ROTATED YET
//on the other hand, the functions of the type _txyz
//only accept a ROTATED xp, so let us not forget about rotating

//question about the ORDER of VELOCITY and PRESSURE
// Here, velocity must go BEFORE pressure
//u_value[0]=ux, u_value[1]=uy, u_value[2]=uz, u_value[3]=up
//therefore, the routine that calls this ic_read (GenIc) has an ORDER in IT

 //this law is given in the NON-DIMENSIONAL domain, with NON-DIMENSIONAL values
 //NONDIMENSIONAL pressure distribution, fundamental!!
 
void EqnNS::ic_read(const double xp[],double u_value[], const double el_xm[]) const {

  const double Uref = _phys.get("Uref");
  double pref = _phys.get("pref");
  const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");

  double udes = _phys.get("udes");
  
  
  Box* box= static_cast<Box*>(_mesh.GetDomain());
  
  double         lb[DIMENSION];
  double         le[DIMENSION];
  lb[0] = box->_lb[0];//already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
#if DIMENSION == 3
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
#endif
  
  
  double x_rotshift[DIMENSION];
  _mesh._domain->TransformPointToRef(xp,x_rotshift); 
//at this point, the coordinates are transformed into the REFERENCE BOX, so you can pass them to the Pressure function

//TODO here you should also rotate the ELEMENT COORDINATES  
  
//rotation of the function  
    double thetaz = box->_domain_rtmap.get("thetaz");

  
#if (DIMENSION==2)

const double magnitude = /*udes**/1.*(x_rotshift[0] - box->_lb[0])*(box->_le[0]-x_rotshift[0] )/( Uref); 
    u_value[0] = -sin(thetaz)*magnitude;
    u_value[1] = cos(thetaz)*magnitude; 

//what if you start... well, after the first state run, you already start with the analytical solution
//well, you should put the DESIRED SOLUTION equal to the ANALYTICAL solution without control, so that u - u_d is IDENTICALLY ZERO!
double press_tmp[1]; 
       press_tmp[0]=0.;
 _eqnmap._qtymap.get_qty("Qty_Pressure")->Function_txyz(0.,x_rotshift,press_tmp);   
 u_value[2] = press_tmp[0];

#elif (DIMENSION==3)

  u_value[0] = 0.;
  u_value[1] = 0.*/*udes**/(x_rotshift[0] - box->_lb[0])*(box->_le[0]-x_rotshift[0] )/( Uref);
  u_value[2] = 0.;

double press_tmp[1]; 
       press_tmp[0]=0.;
  _eqnmap._qtymap.get_qty("Qty_Pressure")->Function_txyz(0.,x_rotshift,press_tmp);
u_value[3]= 0.*press_tmp[0];
 
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  {  u_value[1] = (x_rotshift[0] - lb[0])*(le[0] - x_rotshift[0])*(x_rotshift[2] - lb[2])*(le[2]-x_rotshift[2])/Uref;  }
  
#endif


  return;
}


//=======================
//the implementation of these boundary conditions is related to the particular Domain
//Here we are picking a Box
//So, you get the domain name from the Domain. If it is not a box, you abort.
//the imposition of the boundary conditions is related to the Equation.
//Clearly, it depends on the domain
//So for different domains we would have different parts here, with if's.
//We cannot associate this function to the Box or the Cylinder because 
//it depends on the OPERATORS involved in the EQUATION,
//so it must stay stick to the Equation, which is a bunch of operators
//every application has only one domain, but if you want to use different 
//domains in the same equation you have to specify it here...
//also, changing the domain would mean changing the functions in the Physics User Quantities,
//so in general we do not automatically switch the domain so quickly

//So, for every Domain we have a different implementation 
// of this function
//The idea is: i have to get the Box from where i put it.
//The point is that i set it as a domain but it is also a Box
//So i have to do a CAST from Domain to Box


void EqnNS::elem_bc_read(const double el_xm[],int& surf_id, double value[],int el_flag[]) const {
//el_xm[] is the NON-DIMENSIONAL node coordinate // lb,le are NONDIMENSIONALIZED

const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");

  

Box* box= static_cast<Box*>(_mesh.GetDomain());

  double         lb[DIMENSION];
  double         le[DIMENSION];
  lb[0] = box->_lb[0];//already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
#if DIMENSION == 3
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
#endif
  
  double x_rotshift[DIMENSION];
  _mesh._domain->TransformPointToRef(el_xm,x_rotshift);

 
 #if (DIMENSION==2)
  

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) { //left
surf_id=44; 
     el_flag[NN]=1;
     el_flag[TT]=1;   
  value[NN]=0.;
  value[TT]=0.;/*-4.*/
    
  }


 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0]) -(x_rotshift[0]) < bdry_toll){ //right
surf_id=66; 
     el_flag[NN]=1; 
     el_flag[TT]=1;
       value[NN]=0.;
       value[TT]=0.;/*+4.*/ 
}
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom
surf_id=22; 

      el_flag[NN]=0;    //no normal component
      el_flag[TT]=1;    //yes tangential component
        value[NN]=0.;
        value[TT]=0.;
}
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  { //top
 surf_id=88;

     el_flag[NN]=0;     //no normal component
     el_flag[TT]=1;     //yes tangential component
       value[NN]=0.;
       value[TT]=0.;
    
  }
  
 

#elif (DIMENSION==3)
 
  
  
 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) { //left
surf_id=44;  
     el_flag[NN]=1;    //yes normal component
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) -x_rotshift[0] < bdry_toll){ //right
surf_id=66;
     el_flag[NN]=1;    //yes normal component 
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
   
}
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom

   surf_id=22;

     el_flag[NN]=0;    //no normal component
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }

  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  { //top
surf_id=88;

     el_flag[NN]=0;    //no normal component
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }
  
 if ( x_rotshift[2] > -bdry_toll &&  x_rotshift[2] < bdry_toll ) { //symmetry

surf_id=11;

     el_flag[NN]=1;    //yes normal (equal to zero)
     el_flag[TT]=0;    //no tangential (symmetry)
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }
  if ((le[2]-lb[2]) - x_rotshift[2] > -bdry_toll &&  (le[2]-lb[2]) -x_rotshift[2] < bdry_toll)  {
surf_id=77; 
     el_flag[NN]=1;    //yes normal component //it can be zero also i think
     el_flag[TT]=0;    //no tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the tangential
  value[3]=0.;  //value of the tangential
  
  }


#endif


  return;
}




/// This function  defines the boundary conditions for the NS system:
//TODO ok, here we need the positions of the SCALAR quantities.
//So, since we know that we have one vector quantity and one scalar quantity,
//we must which are called qty velocity and qty pressure,
//we pick their position from themselves
//of course for the vector quantities their positions are are consecutive,
// and also depending on the dimension

//


void EqnNS::bc_read(const double xp[], const double /*normal */[],int bc_flag[]) const {
//xp[] is the NON-DIMENSIONAL node coordinate

  const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");

  //Pick the positions of the scalar quantities
  int pos_ux = _eqnmap._qtymap.get_qty("Qty_Velocity")->_pos;
  int pos_uy = _eqnmap._qtymap.get_qty("Qty_Velocity")->_pos + 1;
  #if DIMENSION == 3
  int pos_uz = _eqnmap._qtymap.get_qty("Qty_Velocity")->_pos + 2;
  #endif
  int pos_up = DIMENSION - 1 + _eqnmap._qtymap.get_qty("Qty_Pressure")->_pos;

  
Box* box = static_cast<Box*>(_mesh.GetDomain());

  double         lb[DIMENSION];
  double         le[DIMENSION];
  lb[0] = box->_lb[0];//already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
#if DIMENSION == 3
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
#endif
  
  double x_rotshift[DIMENSION];
  _mesh._domain->TransformPointToRef(xp,x_rotshift);
  
  
#if (DIMENSION==2)
  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
     bc_flag[pos_ux]=0;
     bc_flag[pos_uy]=0;  
//   bc_flag[pos_up]=0;
  }

 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
    bc_flag[pos_ux]=0;
    bc_flag[pos_uy]=0;
//  bc_flag[pos_up]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
     bc_flag[pos_ux]=0;   //u dot t
//     bc_flag[pos_uy]=0;
    bc_flag[pos_up]=0;
//   if (bc_flag[pos_up] !=0)  bc_flag[2]=1;  //tau dot n   //
                     //remember that this doesnt mean that the pressure at the boundary is fixed!
		     //its not a Dirichlet boundary condition for pressure!
		     //the initial guess of p at the boundary is changed by the equation after one step
		     //none of the pressure nodes are fixed, they are all COMPUTED after the first step
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
    bc_flag[pos_ux]=0;
//     bc_flag[pos_uy]=0;
    bc_flag[pos_up]=0;   //    if (bc_flag[2] !=0)  bc_flag[pos_up]=1;
  }
  

#elif (DIMENSION==3)

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[pos_ux]=0;    //u dot n
    bc_flag[pos_uy]=0;    //u x n
    bc_flag[pos_uz]=0;    //u x n 
//  bc_flag[pos_up]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll ) {  //right of the RefBox
    bc_flag[pos_ux]=0;    //u dot n
    bc_flag[pos_uy]=0;   //u x n
    bc_flag[pos_uz]=0;   //u x n
//  bc_flag[pos_up]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )  {  //bottom  of the RefBox
     bc_flag[pos_ux]=0;      //u x n
//      bc_flag[pos_uy]=0;   //u dot n   //leave this free for VELOCITY INLET
     bc_flag[pos_uz]=0;      //u x n
//      bc_flag[pos_up]=0;     //tau dot n //pressure
  }
  
  if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  {  //top of the  of the RefBox
     bc_flag[pos_ux]=0;     //u x n
//      bc_flag[pos_uy]=0;  //u dot n   //leave this free for outlet
     bc_flag[pos_uz]=0;     //u x n
      bc_flag[pos_up]=0;     //tau dot n  //PRESSURE OUTLET
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
       if (bc_flag[pos_ux] == 1)  bc_flag[pos_ux] = _phys.get("Fake3D");   //u x n  //check it for all equations
       if (bc_flag[pos_uy] == 1)  bc_flag[pos_uy] = _phys.get("Fake3D");   //u x n          //leave this free for 2D
      bc_flag[pos_uz]=0;                                               //u dot n  
//       bc_flag[pos_up]=0;
  }
  
  if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  {
     if (bc_flag[pos_ux] == 1) bc_flag[pos_ux] = _phys.get("Fake3D");      //u x n
     if (bc_flag[pos_uy] == 1) bc_flag[pos_uy] = _phys.get("Fake3D");      //u x n      //leave this free for 2D
     bc_flag[pos_uz] = 0;                                                  //u dot n
//      bc_flag[pos_up]=0;
  }

// // // OPTIMIZATION ANTALYA
// // // I was doing a Dirichlet inlet to have a better control on reverse flow
// // //   // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
// // //   //    boundary conditions box
// // //   if (xp[0] < lxb*ILref + bdry_toll) {
// // //     bc_flag[0]=0;    //u dot n
// // //     bc_flag[1]=0;    //u x n
// // //     bc_flag[2]=0;    //u x n 
// // // //  bc_flag[3]=0;
// // //   }
// // //   if (xp[0] > lxe*ILref - bdry_toll) {
// // //     bc_flag[0]=0;    //u dot n
// // //     bc_flag[1]=0;   //u x n
// // //     bc_flag[2]=0;   //u x n
// // // //  bc_flag[3]=0;
// // //   }
// // //   if (xp[1] < lyb*ILref + bdry_toll) { //  INLET
// // //      bc_flag[0]=0;      //u x n
// // //       bc_flag[1]=0;   //u dot n   //leave this free for inlet
// // //      bc_flag[2]=0;      //u x n
// // // //      bc_flag[3]=0;     //tau dot n //pressure
// // //   }
// // //   if (xp[1] > lye*ILref - bdry_toll) { //  OUTLET
// // //      bc_flag[0]=0;     //u x n
// // // //      bc_flag[1]=0;  //u dot n   //leave this free for outlet
// // //      bc_flag[2]=0;     //u x n
// // //       bc_flag[3]=0;     //tau dot n
// // //   }
// // //   if (xp[2] < lzb*ILref + bdry_toll) {  //current
// // //       if (bc_flag[0] == 1)  bc_flag[0]=_phys.get_par("Fake3D");   //u x n  //check it for all equations
// // //       if (bc_flag[1] == 1)  bc_flag[1]=_phys.get_par("Fake3D");   //u x n          //leave this free for 2D
// // //       bc_flag[2]=0;                                               //u dot n  
// // // //       bc_flag[3]=0;
// // //   }
// // //   if (xp[2] > lze*ILref - bdry_toll) {  //current
// // //      if (bc_flag[0] == 1) bc_flag[0]=_phys.get_par("Fake3D");      //u x n
// // //      if (bc_flag[1] == 1) bc_flag[1]=_phys.get_par("Fake3D");      //u x n      //leave this free for 2D
// // //      bc_flag[2]=0;                                                  //u dot n
// // // //      bc_flag[3]=0;
// // //   }
#endif


   return;
}




///********************MHD*******************
///********************MHD*******************
///********************MHD*******************

//===============================================================
void EqnMHD::ic_read(const double xp[],double u_value[], const double el_xm[]) const {

  const double Uref = _phys.get("Uref");
  const double Bref = _phys.get("Bref");
 
#if (DIMENSION==2)
  u_value[0] = 0./Bref;
  u_value[1] = 0./Bref;
  u_value[2] = 0./(Uref*Bref);

#elif  (DIMENSION==3)
  u_value[0] = 0./Bref;
  u_value[1] = 0./Bref;
  u_value[2] = 0./Bref;
  u_value[3] = 0./(Uref*Bref);
  
#endif

  return;
}


//===============================================================
void EqnMHD::elem_bc_read(const double el_xm[],int& surf_id, double value[],int el_flag[]) const {
//el_xm[] is the NON-DIMENSIONAL node coordinate
// lb,le are NONDIMENSIONALIZED
  //in this way lb,le,el_xm,x_rotshift are ALL nondimensional, so you can compare them!

const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");
      
  

Box* box = static_cast<Box*>(_mesh.GetDomain());

  double         lb[DIMENSION];
  double         le[DIMENSION];
  lb[0] = box->_lb[0];//already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
#if DIMENSION == 3
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
#endif
  
  double x_rotshift[DIMENSION];
  _mesh._domain->TransformPointToRef(el_xm,x_rotshift);

 
 #if (DIMENSION==2)
  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) { //left
surf_id=44;
     el_flag[NN]=1;     //yes normal component //no press integral
     el_flag[TT]=1;     //yes tang component   
  value[NN]=0.;
  value[TT]=0.;

    
  }


 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0]) -(x_rotshift[0]) < bdry_toll){ //right
surf_id=66;
     el_flag[NN]=1; 
     el_flag[TT]=1;
       value[NN]=0.;
       value[TT]=0.;
}
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom
surf_id=22; 

      el_flag[NN]=0;    //no normal component     //b.n free
     el_flag[TT]=1;    //yes tangential component //bxn like uxn
  value[NN]=0.;
  value[TT]=0.;
}
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  { //top
 surf_id=88;

     el_flag[NN]=0;     //b.n free
     el_flag[TT]=1;     //bxn like uxn
  value[NN]=0.;
  value[TT]=0.;
    
  }
  

  

#elif (DIMENSION==3)

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) { //left
surf_id=44;
     el_flag[NN]=1;    //yes normal component
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }

  if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) -x_rotshift[0] < bdry_toll){ //right
surf_id=66;
     el_flag[NN]=1;    //yes normal component 
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
   
}
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom

   surf_id=22;

     el_flag[NN]=0;    //no normal component
     el_flag[TT]=1;    //yes tang component //bxn fixed
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }
//   if (el_xm[1] > lye*ILref - bdry_toll) { //  OUTLET
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  { //top
surf_id=88;

     el_flag[NN]=0;    //no normal component
     el_flag[TT]=1;    //yes tang component  //bxn fixed
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }
  
 if ( x_rotshift[2] > -bdry_toll &&  x_rotshift[2] < bdry_toll ) { //current

surf_id=11;

     el_flag[NN]=1;    //yes normal (equal to zero)(SYMMETRY)
     el_flag[TT]=0;    //no tangential (symmetry)
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }

  if ((le[2]-lb[2]) - x_rotshift[2] > -bdry_toll &&  (le[2]-lb[2]) -x_rotshift[2] < bdry_toll)  {
surf_id=77;
     el_flag[NN]=1;    //yes normal component (SYMMETRY)
     el_flag[TT]=0;    //no tang component(SYMMETRY)
  value[0]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the tangential
  value[3]=0.;  //value of the tangential
  
  }


#endif


  return;
}





/// This function defines the boundary conditions for the MHD system:
    
 //Bxn = (b+Be)xn = 0 would mean perfectly insulating   
 //our present situation is such that:
//on the walls:
//bxn = Bexn = 0 => Bxn = 0
//on the inlet and outlet we cannot speak of WALLS
//how do we set BOUNDARY/INTERFACE conditions on the INLET/OUTLET?
//right out of the boundary there will still be the same material
//so we have to put INFINITY|SYMMETRY conditions
//because we are studying an infinite channel:
//   dB/dy=0
//that would also come out from a formulation with the Laplacian
//So you dont have to fix the test function, because the integrand 
//is specified (equal to zero... but specified)
// So for symmetry you LET FREE all the quantities
// and that's ok. That means that the possible
//BOUNDARY INTEGRALS containing FIRST DERIVATIVES 
//of as a form of NORMAL DERIVATIVE, CURL DERIVATIVE or whatever
//are set TO ZERO, 
//not because the projection onto the normal is ZERO
//but because you assume that ALL THE NORMAL DERIVATIVES ARE ZERO.

//however we have
 //bxn = 0
 //Bexn = either FIXED or CONTROLLED value
void EqnMHD::bc_read(const double xp[], const double /*normal */[],int bc_flag[]) const {

    const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");
  
    
Box* box = static_cast<Box*>(_mesh.GetDomain());

  double         lb[DIMENSION];
  double         le[DIMENSION];
  lb[0] = box->_lb[0];//already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
#if DIMENSION == 3
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
#endif
  
  double x_rotshift[DIMENSION];
  _mesh._domain->TransformPointToRef(xp,x_rotshift);
    
  
#if (DIMENSION==2)

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {
//     bc_flag[0]=0;  //b.n useless with curl curl
     bc_flag[1]=0;  //bxn
//   bc_flag[2]=0;
  }
  
 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){
//    bc_flag[0]=0;  //b.n useless with curl curl
    bc_flag[1]=0;  //bxn
//  bc_flag[2]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  {
     bc_flag[0]=0;  //bxn
//      bc_flag[1]=0;      //leave this free
//          bc_flag[2]=0;     //comment it, don't do the integral
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {
    bc_flag[0]=0;  //bxn
//     bc_flag[1]=0;      //leave this free
//         bc_flag[2]=0;       //comment it, don't do the integral
 }
  
#else

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) { //  INSULATING
    bc_flag[0]=0;   //b.n
    bc_flag[1]=0; //bxn
    bc_flag[2]=0;   //bxn
//  bc_flag[3]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){ //   INSULATING
    bc_flag[0]=0;   //b.n
      bc_flag[1]=0; //bxn
    bc_flag[2]=0;   //bxn
//  bc_flag[3]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {
    bc_flag[0]=0;       //bxn
//     bc_flag[1]=0;    //b.n     //leave this free for inlet
//     bc_flag[2]=0;    //bxn  //WHY ISNT THIS FIXED AS WELL?
//  bc_flag[3]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {
    bc_flag[0]=0;         //bxn
//     bc_flag[1]=0;      //b.n    //leave this free for outlet
//      bc_flag[2]=0;     //bxn   //WHY ISNT THIS FIXED AS WELL?
    //bc_flag[3]=0;
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {   //  CONDUCTING, now insulating
   if (bc_flag[0] == 1) bc_flag[0]=_phys.get("Fake3D");     //bxn
   if (bc_flag[1] == 1) bc_flag[1]=_phys.get("Fake3D");     //bxn      //leave this free for 2D
//      bc_flag[2]=0;                         //b.n 
//     bc_flag[3]=0;
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {   //  CONDUCTING, now insulating
   if (bc_flag[0] == 1) bc_flag[0]=_phys.get("Fake3D");   //bxn
   if (bc_flag[1] == 1) bc_flag[1]=_phys.get("Fake3D");   //bxn     //leave this free for 2D
//  bc_flag[2]=0;                           //b.n 
//     bc_flag[3]=0;
  }
  
#endif


  return;
}



// ===============================

void EqnNSAD::ic_read(const double xp[],double u_value[], const double el_xm[]) const {
//xp[] is the NON-DIMENSIONAL node coordinate  //TODO: reference values

#if (DIMENSION==2)
  u_value[0] = 0.;
  u_value[1] = 0.;
  u_value[2] = 0.;
#else
  u_value[0] = 0.;
  u_value[1] = 0.;
  u_value[2] = 0.;
  u_value[3] = 0.;
#endif

  return;
}



/// This function  defines the boundary conditions for the NS system:

void EqnNSAD::bc_read(const double xp[], const double /*normal */[],int bc_flag[]) const {
//xp[] is the NON-DIMENSIONAL node coordinate

  const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");
  

  
 Box* box = static_cast<Box*>(_mesh.GetDomain());

  double         lb[DIMENSION];
  double         le[DIMENSION];
  lb[0] = box->_lb[0];//already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
#if DIMENSION == 3
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
#endif
  
  double x_rotshift[DIMENSION];
  _mesh._domain->TransformPointToRef(xp,x_rotshift);
  

#if (DIMENSION==2)
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
     bc_flag[0]=0;
     bc_flag[1]=0;
//   bc_flag[2]=0;
  }
  
 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
    bc_flag[0]=0; 
    bc_flag[1]=0;
//  bc_flag[2]=0;
  }

if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
     bc_flag[0]=0;
//   bc_flag[1]=0;   //comment it, because you leave u_y free
//     bc_flag[2]=0;  //adjoint pressure: computing the integral is not correct, because the function p is prescribed
                    //at the boundary, so deltap = 0 at the boundary.
		    //You dont have to consider the symmetry with the direct equation!
                    //so THERE IS NO BOUNDARY INTEGRAL to be computed
		    //  COMMENT bc FOR THE ADJOINT PRESSURE
                    //The problem is that p_old is modified after one nonlinear step,
		    //so p was changing every time!
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
    bc_flag[0]=0;
//  bc_flag[1]=0;  //comment it, because you leave u_y free
//    bc_flag[2]=0;   //  COMMENT bc FOR THE ADJOINT PRESSURE
  }
#else

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
//  bc_flag[3]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
//  bc_flag[3]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
     bc_flag[0]=0;
//      bc_flag[1]=0;    //INSTEAD THIS MUST CORRESPOND TO THE DIRECT
                         //the SPACE of ADJOINT functions is the same as the SPACE for the DIRECT test functions
                         //if you fix this then you dont control...
     bc_flag[2]=0;
//           bc_flag[3]=0;    // COMMENT bc FOR THE ADJOINT PRESSURE
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
     bc_flag[0]=0;
//      bc_flag[1]=0;     //INSTEAD THIS MUST CORRESPOND TO THE DIRECT
     bc_flag[2]=0;
//           bc_flag[3]=0;   // COMMENT bc FOR THE ADJOINT PRESSURE
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
    if (bc_flag[0] == 1) bc_flag[0]=_phys.get("Fake3D");   //u x n
    if (bc_flag[1] == 1) bc_flag[1]=_phys.get("Fake3D");  //u x n             //leave this free for 2D
      bc_flag[2]=0;                                            //u dot n
//      bc_flag[3]=0;
  }
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
    if (bc_flag[0] == 1) bc_flag[0] = _phys.get("Fake3D");   //u x n
    if (bc_flag[1] == 1) bc_flag[1] = _phys.get("Fake3D");  //u x n             //leave this free for 2D
    bc_flag[2]=0;                                                //u dot n
//      bc_flag[3]=0;
  }
#endif


  return;
}

// // ===============================
// TODO must always be in tune with the DIRECT Equation
void EqnNSAD::elem_bc_read(const double el_xm[],int& surf_id, double value[],int el_flag[]) const {
//el_xm[] is the NON-DIMENSIONAL node coordinate

const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");
      
  

Box* box = static_cast<Box*>(_mesh.GetDomain());

  double         lb[DIMENSION];
  double         le[DIMENSION];
  lb[0] = box->_lb[0];//already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
#if DIMENSION == 3
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
#endif
  
  double x_rotshift[DIMENSION];
  _mesh._domain->TransformPointToRef(el_xm,x_rotshift);


 
 #if (DIMENSION==2)
  
  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) { //left
surf_id=44; 
     el_flag[NN]=1;
     el_flag[TT]=1;   
  value[NN]=0.;
  value[TT]=0.;/*-4.*/

  }

 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0]) -(x_rotshift[0]) < bdry_toll){ //right
surf_id=66; 
     el_flag[NN]=1; 
     el_flag[TT]=1;
       value[NN]=0.;
       value[TT]=0.;/*+4.*/ 
}
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom
surf_id=22; 

      el_flag[NN]=0;    //no normal component
      el_flag[TT]=1;    //yes tangential component
        value[NN]=0.;
        value[TT]=0.;
}
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  { //top
 surf_id=88;

     el_flag[NN]=0;     //no normal component
     el_flag[TT]=1;     //yes tangential component
       value[NN]=0.;
       value[TT]=0.;
  }
  
#elif (DIMENSION==3)

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) { //left
surf_id=44;  
     el_flag[NN]=1;    //yes normal component
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) -x_rotshift[0] < bdry_toll){ //right
surf_id=66;
     el_flag[NN]=1;    //yes normal component 
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
   
}
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom

   surf_id=22;

     el_flag[NN]=0;    //no normal component
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }

  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  { //top
surf_id=88;

     el_flag[NN]=0;    //no normal component
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }
  
 if ( x_rotshift[2] > -bdry_toll &&  x_rotshift[2] < bdry_toll ) { //symmetry

surf_id=11;

     el_flag[NN]=1;    //yes normal (equal to zero)
     el_flag[TT]=0;    //no tangential (symmetry)
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }
  if ((le[2]-lb[2]) - x_rotshift[2] > -bdry_toll &&  (le[2]-lb[2]) -x_rotshift[2] < bdry_toll)  {
surf_id=77; 
     el_flag[NN]=1;    //yes normal component //it can be zero also i think
     el_flag[TT]=0;    //no tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the tangential
  value[3]=0.;  //value of the tangential
  
  }


#endif


  return;

}



// // ===============================
// // ============ MHDAD ===================
// // ===============================


/// This function generates the initial conditions for the NS system:

void EqnMHDAD::ic_read(const double xp[],double u_value[], const double el_xm[]) const {
  
#if (DIMENSION==2)
  u_value[0] = 0.;
  u_value[1] = 0.;
  u_value[2] = 0.;
#else
  u_value[0] = 0.;
  u_value[1] = 0.;
  u_value[2] = 0.;
  u_value[3] = 0.;
#endif

  return;
  
}


/// This function  defines the boundary conditions for the NS system:

void EqnMHDAD::bc_read(const double xp[], const double /*normal */[],int bc_flag[]) const {
//xp[] is the NON-DIMENSIONAL node coordinate

  const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");
  
  
Box* box = static_cast<Box*>(_mesh.GetDomain());

  double         lb[DIMENSION];
  double         le[DIMENSION];
  lb[0] = box->_lb[0];//already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
#if DIMENSION == 3
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
#endif
  
  double x_rotshift[DIMENSION];
  _mesh._domain->TransformPointToRef(xp,x_rotshift);  
  

#if (DIMENSION==2)

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
     bc_flag[0]=0;
     bc_flag[1]=0;
//   bc_flag[2]=0;
  }
  
   if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
    bc_flag[0]=0; 
    bc_flag[1]=0;
//  bc_flag[2]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox  
     bc_flag[0]=0;
//      bc_flag[1]=0;      //u dot n leave this free
//      bc_flag[2]=0;      //NO BOUNDARY ADJOINT MHD PRESSURE
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox  
    bc_flag[0]=0;
//     bc_flag[1]=0;      //u dot n leave this free
//     bc_flag[2]=0;      //NO BOUNDARY ADJOINT MHD PRESSURE 
  }
#else

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
//  bc_flag[3]=0;
  }
  
   if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
//  bc_flag[3]=0;
  }
  
     if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox  
     bc_flag[0]=0;
//      bc_flag[1]=0;    //u dot n 
     bc_flag[2]=0;
//           bc_flag[3]=0; //NO BOUNDARY ADJOINT MHD PRESSURE
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox  

     bc_flag[0]=0;
//      bc_flag[1]=0;     //u dot n
     bc_flag[2]=0;
//           bc_flag[3]=0;//NO BOUNDARY ADJOINT MHD PRESSURE
  }
  
if ( x_rotshift[2] > -bdry_toll &&  x_rotshift[2] < bdry_toll ) { //current  
     if (bc_flag[0] == 1) bc_flag[0]=_phys.get("Fake3D");
     if (bc_flag[1] == 1) bc_flag[1]=_phys.get("Fake3D");
//      bc_flag[2]=0;
//      bc_flag[3]=0;
  }
  
  if ((le[2]-lb[2]) - x_rotshift[2] > -bdry_toll &&  (le[2]-lb[2]) -x_rotshift[2] < bdry_toll)  {
    if (bc_flag[0] == 1) bc_flag[0]=_phys.get("Fake3D");
    if (bc_flag[1] == 1) bc_flag[1]=_phys.get("Fake3D");
//    bc_flag[2]=0;
//      bc_flag[3]=0;
  }
#endif


  return;
}


// // ===============================
// TODO must always be in tune with the DIRECT Equation
void EqnMHDAD::elem_bc_read(const double el_xm[],int& surf_id, double value[],int el_flag[]) const{
//el_xm[] is the NON-DIMENSIONAL node coordinate
// lb,le are NONDIMENSIONALIZED
  //in this way lb,le,el_xm,x_rotshift are ALL nondimensional, so you can compare them!

const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");
      

Box* box = static_cast<Box*>(_mesh.GetDomain());

  double         lb[DIMENSION];
  double         le[DIMENSION];
  lb[0] = box->_lb[0];//already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
#if DIMENSION == 3
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
#endif
  
  double x_rotshift[DIMENSION];
  _mesh._domain->TransformPointToRef(el_xm,x_rotshift);

 
 #if (DIMENSION==2)
  
  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) { //left
surf_id=44; 
     el_flag[NN]=1;
     el_flag[TT]=1;   
  value[NN]=0.;
  value[TT]=0.;/*-4.*/

  }

 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0]) -(x_rotshift[0]) < bdry_toll){ //right
surf_id=66; 
     el_flag[NN]=1; 
     el_flag[TT]=1;
       value[NN]=0.;
       value[TT]=0.;/*+4.*/ 
}
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom
surf_id=22; 

      el_flag[NN]=0;    //no normal component
      el_flag[TT]=1;    //yes tangential component
        value[NN]=0.;
        value[TT]=0.;
}
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  { //top
 surf_id=88;

     el_flag[NN]=0;     //no normal component
     el_flag[TT]=1;     //yes tangential component
       value[NN]=0.;
       value[TT]=0.;
  }
  
#elif (DIMENSION==3)

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) { //left
surf_id=44;  
     el_flag[NN]=1;    //yes normal component
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) -x_rotshift[0] < bdry_toll){ //right
surf_id=66;
     el_flag[NN]=1;    //yes normal component 
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
   
}
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom

   surf_id=22;

     el_flag[NN]=0;    //no normal component
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }

  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  { //top
surf_id=88;

     el_flag[NN]=0;    //no normal component
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }
  
 if ( x_rotshift[2] > -bdry_toll &&  x_rotshift[2] < bdry_toll ) { //symmetry

surf_id=11;

     el_flag[NN]=1;    //yes normal (equal to zero)
     el_flag[TT]=0;    //no tangential (symmetry)
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }
  if ((le[2]-lb[2]) - x_rotshift[2] > -bdry_toll &&  (le[2]-lb[2]) -x_rotshift[2] < bdry_toll)  {
surf_id=77; 
     el_flag[NN]=1;    //yes normal component //it can be zero also i think
     el_flag[TT]=0;    //no tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the tangential
  value[3]=0.;  //value of the tangential
  
  }


#endif



  return;

}


///************ MHDCONT_EQUATIONS
///************ MHDCONT_EQUATIONS
///************ MHDCONT_EQUATIONS

void EqnMHDCONT::ic_read(const double xp[], double u_value[], const double el_xm[]) const {

  const double Uref = _phys.get("Uref");
  const double Bref = _phys.get("Bref");
 
#if (DIMENSION==2)
  u_value[0] = Bref/Bref;
  u_value[1] = 0./Bref;
  u_value[2] = 0./(Uref*Bref);

#else
  u_value[0] = Bref/Bref;
  u_value[1] = 0./Bref;
  u_value[2] = 0./Bref;
  u_value[3] = 0./(Uref*Bref);
  
#endif

  return;
}


void EqnMHDCONT::elem_bc_read(const double el_xm[],int& surf_id, double value[],int el_flag[]) const {
//el_xm[] is the NON-DIMENSIONAL node coordinate
// lb,le are NONDIMENSIONALIZED
  //in this way lb,le,el_xm,x_rotshift are ALL nondimensional, so you can compare them!

const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");
      
  

Box* box= static_cast<Box*>(_mesh.GetDomain());

  double         lb[DIMENSION];
  double         le[DIMENSION];
  lb[0] = box->_lb[0];//already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
#if DIMENSION == 3
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
#endif
  
  double x_rotshift[DIMENSION];
  _mesh._domain->TransformPointToRef(el_xm,x_rotshift);
  
 
 #if (DIMENSION==2)
  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) { //left
surf_id=44;
     el_flag[NN]=1;     //yes normal component //no press integral
     el_flag[TT]=1;     //yes tang component   
  value[NN]=0.;
  value[TT]=0.;

    
  }


 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0]) -(x_rotshift[0]) < bdry_toll){ //right
surf_id=66;
     el_flag[NN]=1; 
     el_flag[TT]=1;
       value[NN]=0.;
       value[TT]=0.;
}
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom
surf_id=22; 

      el_flag[NN]=0;    //no normal component     //b.n free
     el_flag[TT]=1;    //yes tangential component //bxn like uxn
  value[NN]=0.;
  value[TT]=0.;
}
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  { //top
 surf_id=88;

     el_flag[NN]=0;     //b.n free
     el_flag[TT]=1;     //bxn like uxn
  value[NN]=0.;
  value[TT]=0.;
    
  }
  

  

#elif (DIMENSION==3)

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) { //left
surf_id=44;
     el_flag[NN]=1;    //yes normal component
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }

  if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) -x_rotshift[0] < bdry_toll){ //right
surf_id=66;
     el_flag[NN]=1;    //yes normal component 
     el_flag[TT]=1;    //yes tang component
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
   
}
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom

   surf_id=22;

     el_flag[NN]=0;    //no normal component
     el_flag[TT]=1;    //yes tang component //bxn fixed
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }

  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  { //top
surf_id=88;

     el_flag[NN]=0;    //no normal component
     el_flag[TT]=1;    //yes tang component  //bxn fixed
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }
  
 if ( x_rotshift[2] > -bdry_toll &&  x_rotshift[2] < bdry_toll ) { //current

surf_id=11;

     el_flag[NN]=1;    //yes normal (equal to zero)(SYMMETRY)
     el_flag[TT]=0;    //no tangential (symmetry)
  value[NN]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the  tangential
  value[3]=0.;  //value of the tangential
  }

  if ((le[2]-lb[2]) - x_rotshift[2] > -bdry_toll &&  (le[2]-lb[2]) -x_rotshift[2] < bdry_toll)  {
surf_id=77;
     el_flag[NN]=1;    //yes normal component (SYMMETRY)
     el_flag[TT]=0;    //no tang component(SYMMETRY)
  value[0]=0.;  //value of the normal 
  value[1]=0.;  //value of the tangential
  value[2]=0.;  //value of the tangential
  value[3]=0.;  //value of the tangential
  
  }


#endif



  return;
}




/// This function  defines the boundary conditions for the MHD system:


void EqnMHDCONT::bc_read(const double xp[], const double /*normal */[],int bc_flag[]) const {

  const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");

  
Box* box= static_cast<Box*>(_mesh.GetDomain());

  double         lb[DIMENSION];
  double         le[DIMENSION];
  lb[0] = box->_lb[0];
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
#if DIMENSION == 3
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
#endif
  
  double x_rotshift[DIMENSION];
  _mesh._domain->TransformPointToRef(xp,x_rotshift);  
  

 
#if (DIMENSION==2)
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) { //left
     bc_flag[0]=0;
     bc_flag[1]=0;
//   bc_flag[2]=0;
  }
  
  if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){  //right
    bc_flag[0]=0; 
    bc_flag[1]=0;
//  bc_flag[2]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
     bc_flag[0]=0;
     bc_flag[1]=0;  //here when you control you must fix it,even if the corresponding b is let free, otherwise the control may use this as well
                      //well,wait,it depends: if you are using the Laplacian for MHD state, then Becont bc's must be consistent with the Laplacian
		      // given there.
		      //if you are using curlxcurl, then Becont should be consistent with curl-curl
		      //on the other hand, Be has a gamma*Laplacian FOR HERSELF... !
		      //so the Becont BC's should be consistent BOTH WITH MHD EQUATION and with MHDCONT EQUATION!
		      
		      
// wait, if you fix all dirichlet, then things go into Becontp pressure
//instead if you fix one Dir and one Neu then nothing goes into pressure

//          bc_flag[2]=0;     //comment it, don't do the integral
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {//top of the  of the RefBox
    
     bc_flag[0]=_phys._physrtmap.get("UseControl");  ////////////////
     bc_flag[1]=_phys._physrtmap.get("UseControl");  ///////////////
//         bc_flag[2]=0;       //comment it, don't do the integral

 }
  
#else

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
//  bc_flag[3]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
//  bc_flag[3]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
//  bc_flag[3]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox
    bc_flag[0]=_phys.get("UseControl");    // 
    bc_flag[1]=_phys.get("UseControl");   //
    
    bc_flag[2]=0;
    //bc_flag[3]=0;
    
  }
   if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
//     bc_flag[3]=0;
  }
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
//     bc_flag[3]=0;
  }
  
#endif //DIMENSION



  return;
}


} //end namespace femus


