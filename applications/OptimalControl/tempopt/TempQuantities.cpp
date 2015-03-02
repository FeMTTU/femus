
//C++ includes
#include <cmath>
#include <iostream>

//library includes
#include "Typedefs.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "Box.hpp"

//application
#include "OptLoop.hpp"
#include "TempQuantities.hpp"

//=================== BEGIN CONSTRUCTORS ================================
// ==================================================================
// ==================================================================


//===========================================================================
Temperature::Temperature(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in) 
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { }

//===========================================================================
TempLift::TempLift(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in, const TimeLoop & timeloop_in) : _my_timeloop(timeloop_in), 
 Quantity(name_in,qtymap_in,dim_in,FEord_in) { }

//===========================================================================
TempAdj::TempAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in) 
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { }

//===========================================================================
TempDes::TempDes(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in) 
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { }

//========================
Pressure::Pressure(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {
  
  for (uint i=0;i<dim_in;i++) _refvalue[i] = qtymap_in.GetInputParser()->get("pref");
}


//=========================================================================
Velocity::Velocity(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {  

   for (uint i=0;i<dim_in;i++) _refvalue[i] =  qtymap_in.GetInputParser()->get("Uref");
}


//========================
Pressure2::Pressure2(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

  for (uint i=0;i<dim_in;i++) _refvalue[i] = qtymap_in.GetInputParser()->get("pref");
}

//=================== END CONSTRUCTORS ================================
// ==================================================================
// ==================================================================


//=============================================================
///analytical velocity for Hartmann flow
// difference between get_par and optsys:
// in both cases you are "dynamic" somehow

void Velocity::Function_txyz(const double /*t*/,const double* xp, double* func) const {

  Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());
  // we should do this static_cast in the QUANTITY or QUANTITY MAP constructor
  //if there is some domain shape, we see what type it is and we do the static cast
  //if there is no domain shape, we dont need the domain.
  
    //=====ROTATION of the Function
    //this one is about the reference frame, 
    //here instead we dont want to change the reference frame, 
    //we just want to rotate the function but in a straight reference
  const double thetaz = 0.*3.14/2.;  //const double thetaz = box->_domain_rtmap.get("thetaz");
  
  //====== Physics
//   TempPhysics *optphys; optphys = static_cast<TempPhysics*>(&(_qtymap._phys));
  
  const double rhof   = _qtymap.GetInputParser()->get("rho0");
  const double Uref   = _qtymap.GetInputParser()->get("Uref");
  const double Lref   = _qtymap.GetInputParser()->get("Lref");

  const double DpDz   = 1./*0.5*/;  //AAA: change it according to the pressure distribution!!!

  double DpDzad = DpDz*Lref/(rhof*Uref*Uref);

//   double Re  = optphys->_Re;


//   double Lhalf = 0.5*(box->_le[0] - box->_lb[0]);
//   double Lmid  = 0.5*(box->_le[0] + box->_lb[0]);

//   double xtr = xp[0] - Lmid/*/Lref*/;

  
 //constant for the real reference length in the Hartmann number
//   const double LHm =2.;   //this is because the reference length for Hm is HALF THE WIDTH of the domain, which is Lref=1 now
  const double magnitude = 5.*DpDzad*xp[0]*(1. - xp[0]);
 
  
  func[0] = -sin(thetaz)*magnitude/*/Uref*/;
  func[1] = cos(thetaz)*magnitude;
                                       //add a 4. to the denominator
				       //should check the difference between L and Lref
                                       //TODO check this nondimensionalization
    if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  func[2] = 0./*/Uref*/;
    }


  return;

 
}



//============================================================= 
void Velocity::strain_txyz(const double /*t*/, const double* xyz,double strain[][DIMENSION]) const {

//here, tau is a tensor, so tau dot n is a vector which in general has a NORMAL and a TANGENTIAL component  
  
    const double Lref = _qtymap.GetInputParser()->get("Lref");
      double ILref = 1./Lref;
      const double lye = _qtymap.GetMeshTwo()->GetDomain()->_domain_rtmap.get("lye");
//   const double x=xyz[0];
  const double y=xyz[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  const double z=xyz[2];
  }
  
  strain[0][0] = 0.;                     //ux,x
  strain[0][1] = strain[1][0] = 0. ;//0.5*(uy,x+ux,y) 
  strain[1][1] = 0.*(-(lye*ILref-y));                    //uy,y
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  strain[0][2] = strain[2][0] = 0. ;  //0.5*(uz,x+ux,z) 
  strain[1][2] = strain[2][1] = 0. ;  //0.5*(uy,z+uz,y)                                     
  strain[2][2] = 0. ;                    //uz,z
  }

return;
}


//=============================================================
/// prescribed pressure at the boundary
//no initial condition for pressure is required, because it has no time derivative
//only the boundary condition in the Neumann part of the boundary has to be enforced
//so there is also no problem about consistency between IC and BC values,
// because there are NO IC values for pressure! 

void Pressure::Function_txyz(const double t, const double* xp,double* func) const {
  //the coordinates (x,y,z,t) of the VOLUME domain are NON-dimensional
  //and the function value must be nondimensional as well
  //this function receives an ALREADY ROTATED COORDINATE!

  ///this function may receive the values at the DOFS or at the GAUSS POINTS also;
  ///in both cases, the coordinates must be in the REFBOX DOMAIN,
  ///so use xyzrefbox._val_g or xyzrefbox._val_dofs


//here you see that you perform a casting
//of BOTH the DOMAIN and the PHYSICS
//this is clear because the Base functions
//only yield the Base datatypes,
//so if you are in a CHILD CLASS you have to 
//do the CASTS towards the CHILD classes.
//Clearly, the goal would be to do the castings
//NOT INSIDE SMALL ROUTINES but as class members 
//or something, anyway at higher level,
//so that these small functions do not need to do it repeatedly

//the point is that while one of the two casts can be done "a priori"
//because one may establish "a priori" the set of specific Domain Shapes
//Box,Cylinder, whatever  (well, one could actually add it, and update the library...)
//the cast of the OptPhys cannot be done automatically
// so it should be done
//    in ALL the SPECIFIC Quantity constructors
//and in ALL the SPECIFIC Equation constructors
//since all these classes are application - specific,
//you do not perturb the library.

// also, we have to think how to do with the Domains, because we do not want 
// the user to need to update the library for every different domain
// we must think of an Application Specific domain 
// that must be cast after its introduction in the Mesh class.
// Everyone can reach the Domain through the mesh class as a FATHER domain;
//how can we convert it to a specific domain EXPLICITLY
// and STILL STAYING in the Mesh class WITHOUT PERTURBING it?
//the basic classes Mesh, MultiLevelProblemTwo, QuantityMap handle FATHER THINGS.
//the application-specific classes (Equation, Quantity, Physics) handle CHILD things.
//in this passage we must CONVERT at the highest possible level.
//we should INSTANTIATE the CHILDREN in the applications 
// and PASS THEM SEPARATELY to the basic classes and application-specific classes
// no we cant do like this, we pass only once to the basics and then
//convert in the app specific ones.

//yes, for the domain shapes we must find a way to avoid the switch Box Cylinder,
// because if one adds a new shape one SHOULD MODIFY the GENCASE also 
//and we'd want to keep it as small as possible.
// The point is that for now the Gencase CANNOT LIVE without SPECIFIC BOX information
// because we implemented the functions for BOX GENERATION.
//So if one has cylinder one should add the cylinder generation functions
// if one has a Xmas tree one should add the christmas tree generation functions...
// NO, we clearly cannot do that
// i would want the gencase to be AS GENERAL as POSSIBLE

Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());  //already nondimensionalized
 
  func[0] =  1./ _qtymap.GetInputParser()->get("pref")*( (box->_le[1] - box->_lb[1]) - xp[1] )*(cos(6.*0.*t));

//this equation is in the reference frame CENTERED AT (0,0,0)  
  
  return;
  }


  
  
//===============
void TempDes::Function_txyz(const double /*t*/, const double* /*xp*/,double* temp) const {
  
  temp[0] = 0.9;
 
  return;
}


//===============
void TempAdj::Function_txyz(const double/* t*/, const double* /*xp*/,double* temp) const {
  
  temp[0] = 0.;
 
  return;
}
  
// =================================================
void TempLift::Function_txyz(const double /*t*/, const double* xp,double* temp) const {

  const double Tref = _qtymap.GetInputParser()->get("Tref");

  Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());  
  
  temp[0] = 100.*(xp[0])*( ( box->_le[0] - box->_lb[0]) - xp[0])*(xp[1])*( ( box->_le[1] - box->_lb[1]) - xp[1])/Tref;
 
  
  return;
  }

// =================================================
void Temperature::Function_txyz(const double/* t*/, const double* xp,double* temp) const {

  const double Tref = _qtymap.GetInputParser()->get("Tref");

  Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());  
  
  temp[0] = 100.*(xp[0])*( ( box->_le[0] - box->_lb[0]) - xp[0])/Tref;
 
  
  return;
  }
  
  
// =================================================
  //the coordinates (x,y,z,t) of the VOLUME domain are NON-dimensional
  //and the function value must be nondimensional as well
 //-----Nonhomogeneous Neumann-------
 // Qflux = - k grad(T) by definition
//  QfluxDOTn>0: energy flows outside (cooling)  QfluxDOTn<0: energy flows inside (heating)
void Temperature::heatflux_txyz(const double /*t*/, const double* /*xyz*/, double* qflux) const {

// std::cout << "Temperature: Heatflux, check which coordinates are passed in here" << std::endl;
//     Box* box= static_cast<Box*>(_qtymap._phys._mesh->GetDomain());
//   const double thetaz = box->_domain_rtmap.get("thetaz");

     qflux[0]=-2.1*0./**cos(thetaz)*/;
     qflux[1]=0./**sin(thetaz)*/;
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
      qflux[2]=0.;
    }

  return;
  }

// =================================================
void Pressure2::Function_txyz(const double/* t*/, const double* xp,double* value) const {

  value[0] = 1.;
  
  return;
  
  }
  
  

void Velocity::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  
Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);

if (_qtymap.GetMeshTwo()->get_dim() == 2) {
  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
     bc_flag[0]=0;
     bc_flag[1]=0;  
  }

 if ( (le[0]-lb[0]) - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll  ){ //right of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
     bc_flag[0]=0;
     bc_flag[1]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox

//outflow only on part of the outlet
  if ( (x_rotshift[0]) > 0.70*(le[0]-lb[0]) ) {
      bc_flag[0]=0;      //ux
//     bc_flag[1]=0;     //uy
 }  //end part outflow
    else {  
      bc_flag[0]=0;  //ux
      bc_flag[1]=0;   //uy
    }
  
  } //top RefBox
  
}  //end dim 2

  else if (_qtymap.GetMeshTwo()->get_dim() == 3) {

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;    //u dot n
    bc_flag[1]=0;    //u x n
    bc_flag[2]=0;    //u x n 
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox
    bc_flag[0]=0;    //u dot n
    bc_flag[1]=0;   //u x n
    bc_flag[2]=0;   //u x n
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
     bc_flag[0]=0;      //u x n
//      bc_flag[1]=0;   //u dot n   //leave this free for inlet
     bc_flag[2]=0;      //u x n
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
     bc_flag[0]=0;     //u x n
//      bc_flag[1]=0;  //u dot n   //leave this free for outlet
     bc_flag[2]=0;     //u x n
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
//       if (bc_flag[0] == 1)  bc_flag[0]=_phys.get_par("Fake3D");   //u x n  //check it for all equations
//       if (bc_flag[1] == 1)  bc_flag[1]=_phys.get_par("Fake3D");   //u x n          //leave this free for 2D
      bc_flag[2]=0;                                               //u dot n  
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
//      if (bc_flag[0] == 1) bc_flag[0]=_phys.get_par("Fake3D");      //u x n
//      if (bc_flag[1] == 1) bc_flag[1]=_phys.get_par("Fake3D");      //u x n      //leave this free for 2D
      bc_flag[2]=0;                                                  //u dot n
  }
  
}  //end dim 3

  return;
 
}

  
  
void Pressure::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  
Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);


  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
//      bc_flag[0]=0;
  }

 if ( (le[0]-lb[0]) - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll  ){ //right of the RefBox
//     bc_flag[0]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
//      bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox

//outflow only on part of the outlet
  if ( (x_rotshift[0]) > 0.70*(le[0]-lb[0]) ) {
      bc_flag[0]=0;
 }  //end part outflow
    else {  
//       bc_flag[0]=0;
    }
  
  } //top RefBox
  

  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
//      bc_flag[0]=0;
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
//      bc_flag[0]=0;
  }

  } // end dim 3

  return;
 
}  
  
  
 void Temperature::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
// T' and its adjoint must be Dirichlet homogeneous everywhere on the boundary, by definition.


  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  

  Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);
  
  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {  //left of the RefBox
      bc_flag[0]=0; //always fixed
  }

  if ( (le[0]-lb[0])  - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll)  { //right of the RefBox
      bc_flag[0]=0; //always fixed
   }
   
  if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
      bc_flag[0]=0; //always fixed
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox
      bc_flag[0]=0; //always fixed
  }


    if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
     bc_flag[0]=0;
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
      bc_flag[0]=0;
  }
  
 } //end dim 3
  

  
  return;
}

 void TempLift::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
// T' and its adjoint must be Dirichlet homogeneous everywhere on the boundary, by definition.


  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  

  Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);
  

  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {  //left of the RefBox
  
    if  (_my_timeloop._curr_t_idx < 1)    { bc_flag[0] = 0; }
    else if ( (x_rotshift[1]) < 0.4*(le[1]-lb[1])  ||  (x_rotshift[1]) > 0.6*(le[1]-lb[1]) )  {  bc_flag[0]=0;  } 
 
  }

 if ( (le[0]-lb[0])  - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll)  { //right of the RefBox
     bc_flag[0]=0;
   }

  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox

     //===== START FROM A SIMPLE STATE SOLUTION =========
     if  (_my_timeloop._curr_t_idx < 1)   { bc_flag[0] = 0; }  //=====CONTROL//////////////
      //===== START FROM A SIMPLE STATE SOLUTION =========
     else if  ( (x_rotshift[0]) < 0.25*(le[0] - lb[0]) || ( x_rotshift[0]) > 0.75*(le[0] - lb[0]) )  {  bc_flag[0]=0; }
     }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox
      bc_flag[0]=0;
    }


  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
     bc_flag[0]=0;
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
      bc_flag[0]=0;
  }
  
  } //end dim 3
  
  
  return;
}


 void TempAdj::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {

  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  

  Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);


  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;  //always fixed
  }

  if ( (le[0]-lb[0])  - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll)  { //right of the RefBox
   bc_flag[0]=0; //always fixed
   }
  if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
   bc_flag[0]=0; //always fixed
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox
      bc_flag[0]=0; //always fixed
  }
  

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;
  }
  
  if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox
    bc_flag[0]=0;
  }
  
  if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
     bc_flag[0]=0;    
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
     bc_flag[0]=0;
  }
  
  
  
  
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
     bc_flag[0]=0;
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
      bc_flag[0]=0;
  }
  
  } //end dim 3
  
  return;
}


 void Pressure2::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {

  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  

  Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);


  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {  //left of the RefBox
   bc_flag[0]=0;
  }

  
  return;
}


// ====================== INITIALIZE functions ===================================

void Temperature::initialize_xyz(const double* xp, std::vector< double >& value) const {

      value[0] = 0.; 
  
  return;
}

void TempAdj::initialize_xyz(const double* xp, std::vector< double >& value) const {

      value[0] = 0.; 
      
  return;
}

void TempLift::initialize_xyz(const double* xp, std::vector< double >& value) const {
  
  
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");

  Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]); 

    value[0] = 1.;
  if ((box->_le[1]-box->_lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (box->_le[1]-box->_lb[1]) -(x_rotshift[1]) < bdry_toll)  {
    value[0] = 0.; 
  }


  return;
}

void Velocity::initialize_xyz(const double* xp, std::vector< double >& value) const {
  
  //====== Physics
   const double Uref = _qtymap.GetInputParser()->get("Uref");
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");

 
  Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]); 
//at this point, the coordinates are transformed into the REFERENCE BOX, so you can pass them to the Pressure function

//rotation of the function  
  const double thetaz = box->_domain_rtmap.get("thetaz");

  const double magnitude = 0. /*1.5*(x_rotshift[0] - box->_lb[0])*(box->_le[0]-x_rotshift[0] )/( Uref)*/;
    value[0] = -sin(thetaz)*magnitude;
    value[1] = cos(thetaz)*magnitude; 

//    if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
//    u_value[0] =0;  u_value[1] = 0; }
//  if ( (box->_le[0] - box->_lb[0])  - (x_rotshift[0]) > -bdry_toll && (box->_le[0] - box->_lb[0])  -(x_rotshift[0]) < bdry_toll)  { //right of the RefBox
//       u_value[0] =0;  u_value[1] = 0; }
   
    //==================================
    if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox

//below, inlet
  if ( (x_rotshift[0]) < 0.25*(box->_le[0] - box->_lb[0]) || ( x_rotshift[0]) > 0.75*(box->_le[0] - box->_lb[0]) ) { 
    value[0] = 0.;
    value[1] = 0.;
  }
  else {
    value[0] = 0.; 
    value[1] = 1.; 
    }
  
  }
//============================================

//========================================
//left, inlet
 if  ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {
 
 if ( (x_rotshift[1]) > 0.4*(box->_le[1] - box->_lb[1]) && ( x_rotshift[1]) < 0.6*(box->_le[1]-box->_lb[1]) )  {  //left of the refbox
       value[0] = _qtymap.GetInputParser()->get("injsuc");    
       value[1] = 0.; 
      }
   }   
//============================================

//============================================
//====== outlet
  if ((box->_le[1]-box->_lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (box->_le[1]-box->_lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox

 if ( (x_rotshift[0]) < 0.71*(box->_le[0] - box->_lb[0]) ) {
   value[0] = 0.;
   value[1] = 0.; 
}

  }
  
  
  
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  value[2] = 0.;
  }



  return;
}



void Pressure::initialize_xyz(const double* xp, std::vector< double >& value) const {
  
  Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]); 
//at this point, the coordinates are transformed into the REFERENCE BOX, so you can pass them to the Pressure function

    Function_txyz(0.,&x_rotshift[0],&value[0]);   
  
  return;
}

void Pressure2::initialize_xyz(const double* xp, std::vector< double >& value) const {

  value[0] =  72.*(xp[0]);
  
  return;
}
