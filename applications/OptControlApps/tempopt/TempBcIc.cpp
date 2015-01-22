//C++
#include <cmath>

//library headers
#include "Box.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "MultiLevelMeshTwo.hpp" 
#include "NormTangEnum.hpp"
#include "TimeLoop.hpp"

//application headers
#include "Temp_conf.hpp"
#include "TempQuantities.hpp"
#include "EqnNS.hpp"
#include "EqnT.hpp"


void EqnT::ic_read(const double xp[],double u_value[],const double el_xm[]) const {

//     const double Tref = _phys._physrtmap.get("Tref");
  const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");


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



void EqnT::bc_read(const double xp[],const double /*normal */[],int bc_flag[]) const {
// T' and its adjoint must be Dirichlet homogeneous everywhere on the boundary, by definition.


  const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");
  

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
  
    if  (_my_timeloop._curr_t_idx < 1)  
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
if  (_my_timeloop._curr_t_idx < 1)   bc_flag[1]=0;     //      bc_flag[1]=0; //=====CONTROL//////////////
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



/// This function generates the initial conditions for the NS system:
//xp[] is the NON-dimensional node coordinate
//when this function is called,
//the domain has been NONdimensionalized
//but NOT ROTATED YET
//on the other hand, the functions of the type _txyz
//only accept a ROTATED xp, so let us not forget about rotating

//question about the ORDER of VELOCITY and PRESSURE
// Here, velocity must go BEFORE pressure
//u_value[0]=ux, u_value[1]=uy, u_value[2]=uz, u_value[3]=up
//therefore, the routine that calls this ic_read (GenIc) has an ORDER in IT

 //this law is given in the NON-dimensional domain, with NON-dimensional values
 //NONdimensional pressure distribution, fundamental!!
 
void EqnNS::ic_read(const double xp[],double u_value[],const double el_xm[]) const {

  //====== Physics
   const double Uref = _phys.get("Uref");
  const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");

 
  Box* box= static_cast<Box*>(_mesh.GetDomain());
  
  double* x_rotshift = new double[_mesh.get_dim()];
  _mesh._domain->TransformPointToRef(xp,x_rotshift); 
//at this point, the coordinates are transformed into the REFERENCE BOX, so you can pass them to the Pressure function

//rotation of the function  
  const double thetaz = box->_domain_rtmap.get("thetaz");

  
#if (DIMENSION==2)

const double magnitude = 0. /*1.5*(x_rotshift[0] - box->_lb[0])*(box->_le[0]-x_rotshift[0] )/( Uref)*/;
    u_value[0] = -sin(thetaz)*magnitude;
    u_value[1] = cos(thetaz)*magnitude; 

//    if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
//    u_value[0] =0;  u_value[1] = 0; }
//  if ( (box->_le[0] - box->_lb[0])  - (x_rotshift[0]) > -bdry_toll && (box->_le[0] - box->_lb[0])  -(x_rotshift[0]) < bdry_toll)  { //right of the RefBox
//       u_value[0] =0;  u_value[1] = 0; }
   
    //==================================
    if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox

//below, inlet
  if ( (x_rotshift[0]) < 0.25*(box->_le[0] - box->_lb[0]) || ( x_rotshift[0]) > 0.75*(box->_le[0] - box->_lb[0]) ) { u_value[0] =0;  u_value[1] = 0;}
  else {u_value[0] =0;  u_value[1] = 1.; }
     }
//============================================

//========================================
//left, inlet
 if  ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {
 
 if ( (x_rotshift[1]) > 0.4*(box->_le[1] - box->_lb[1]) && ( x_rotshift[1]) < 0.6*(box->_le[1]-box->_lb[1]) )  {  //left of the refbox
       u_value[0] = _phys.get("injsuc");    u_value[1] = 0; 
      }
   }   
//============================================

//============================================
//====== outlet
  if ((box->_le[1]-box->_lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (box->_le[1]-box->_lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox

 if ( (x_rotshift[0]) < 0.71*(box->_le[0] - box->_lb[0]) ) { u_value[0] = 0.;    u_value[1] = 0; }

  }

//============================================
   
//what if you start... well, after the first state run, you already start with the analytical solution
//well, you should put the DESIRED SOLUTION equal to the ANALYTICAL solution without control, so that u - u_d is IDENTICALLY ZERO!
double press_tmp[1]; 
       press_tmp[0]=0.;
 _eqnmap._qtymap.get_qty("Qty_Pressure")->Function_txyz(0.,x_rotshift,press_tmp);   
 u_value[2] = press_tmp[0];

#elif (DIMENSION==3)

  u_value[0] = 0.;
  u_value[1] = 0.*(x_rotshift[0] - box->_lb[0])*(box->_le[0]-x_rotshift[0] )/( Uref);
  u_value[2] = 0.;

double press_tmp[1]; 
       press_tmp[0] = 0.;
  _eqnmap._qtymap.get_qty("Qty_Pressure")->Function_txyz(0.,x_rotshift,press_tmp);
u_value[3]=press_tmp[0];
  
// if (xp[1] < lyb + bdry_toll) {  u_value[1] = (xp[0] - lxb)*(lxe-xp[0])*(xp[2] - lzb)*(lze-xp[2])/Uref;

  
  #endif

  delete[] x_rotshift;

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
//el_xm[] is the NON-dimensional node coordinate // lb,le are NONdimensionalized

const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");

  

Box* box= static_cast<Box*>(_mesh.GetDomain());

  double* lb = new double[_mesh.get_dim()];
  double* le = new double[_mesh.get_dim()];
  lb[0] = box->_lb[0];//already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_mesh.get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  double* x_rotshift = new double[_mesh.get_dim()];
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


  delete[] lb;
  delete[] le;
  delete[] x_rotshift;

 return;
}




/// This function  defines the boundary conditions for the NS system:

void EqnNS::bc_read(const double xp[],const double /*normal */[],int bc_flag[]) const {
//xp[] is the NON-dimensional node coordinate

  const double bdry_toll = _mesh.GetRuntimeMap().get("bdry_toll");
  
  
Box* box = static_cast<Box*>(_mesh.GetDomain());

  double* lb = new double[_mesh.get_dim()];
  double* le = new double[_mesh.get_dim()];
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_mesh.get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  double* x_rotshift = new double[_mesh.get_dim()];
  _mesh._domain->TransformPointToRef(xp,x_rotshift);
  
  

#if (DIMENSION==2)
  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
     bc_flag[0]=0;
     bc_flag[1]=0;  
//   bc_flag[2]=0;
  }

 if ( (le[0]-lb[0]) - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll  ){ //right of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
//  bc_flag[2]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox

//   if ( (x_rotshift[0]) > 0.25*(le[0]-lb[0]) && ( x_rotshift[0]) < 0.75*(le[0]-lb[0]) ) {
     bc_flag[0]=0;   //u dot t
     bc_flag[1]=0;
//   }
//  bc_flag[2]=0;
//   if (bc_flag[2] !=0)  bc_flag[2]=1;  //tau dot n   //
                     //remember that this doesnt mean that the pressure at the boundary is fixed!
		     //its not a Dirichlet boundary condition for pressure!
		     //the initial guess of p at the boundary is changed by the equation after one step
		     //none of the pressure nodes are fixed, they are all COMPUTED after the first step
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox

//outflow only on part of the outlet
  if ( (x_rotshift[0]) > 0.70*(le[0]-lb[0]) ) {

  
      bc_flag[0]=0;      //ux
//     bc_flag[1]=0;     //uy
      bc_flag[2]=0;      //pressure
//    if (bc_flag[2] !=0)  bc_flag[2]=1;



 }  //end part outflow
    else {  

      bc_flag[0]=0;  //ux
      bc_flag[1]=0;   //uy
//       bc_flag[2]=0;   //pressure

    }
  
  } //top RefBox
  

#elif (DIMENSION==3)

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;    //u dot n
    bc_flag[1]=0;    //u x n
    bc_flag[2]=0;    //u x n 
//  bc_flag[3]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox
    bc_flag[0]=0;    //u dot n
    bc_flag[1]=0;   //u x n
    bc_flag[2]=0;   //u x n
//  bc_flag[3]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
     bc_flag[0]=0;      //u x n
//      bc_flag[1]=0;   //u dot n   //leave this free for inlet
     bc_flag[2]=0;      //u x n
     bc_flag[3]=0;     //tau dot n //pressure
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
     bc_flag[0]=0;     //u x n
//      bc_flag[1]=0;  //u dot n   //leave this free for outlet
     bc_flag[2]=0;     //u x n
      bc_flag[3]=0;     //tau dot n
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
//       if (bc_flag[0] == 1)  bc_flag[0]=_phys.get_par("Fake3D");   //u x n  //check it for all equations
//       if (bc_flag[1] == 1)  bc_flag[1]=_phys.get_par("Fake3D");   //u x n          //leave this free for 2D
      bc_flag[2]=0;                                               //u dot n  
//       bc_flag[3]=0;
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
//      if (bc_flag[0] == 1) bc_flag[0]=_phys.get_par("Fake3D");      //u x n
//      if (bc_flag[1] == 1) bc_flag[1]=_phys.get_par("Fake3D");      //u x n      //leave this free for 2D
     bc_flag[2]=0;                                                  //u dot n
//      bc_flag[3]=0;
  }

#endif

  delete[] lb;
  delete[] le;
  delete[] x_rotshift;
  
  return;
}




//E' chiaro che le condizioni al contorno dipendono dalla geometria e da come implemento l'equazione.
//L'implementazione e' esplicita e viene fatta a priori, quindi sai gia' a priori di quanti blocchi 
// e' composta la matrice di elemento... quindi dentro l'equazione non puoi cambiare il numero 
// di grandezze incognite, A MENO CHE NON FAI DEGLI IFDEF. Piuttosto puoi cambiare quale equazione usi,
// in teoria e' piu' pulito
