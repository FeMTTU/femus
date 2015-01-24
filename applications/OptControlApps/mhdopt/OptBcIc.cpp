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



///COMMENTS on the BC for MHD ====================

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




} //end namespace femus


