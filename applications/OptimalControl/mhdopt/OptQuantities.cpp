
//C++ includes
#include <cmath>
#include <iostream>

//library includes
#include "Typedefs.hpp"
#include "MultiLevelMeshTwo.hpp"

#include "Box.hpp"

//application
#include "OptLoop.hpp"
#include "OptQuantities.hpp"


//=================== BEGIN CONSTRUCTORS ================================
// ==================================================================
// ==================================================================

//===========================================================================
MagnFieldHomX::MagnFieldHomX(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

 for (uint i=0;i<dim_in;i++) _refvalue[i]= qtymap_in.GetInputParser()->get("Bref");
}

//===========================================================================
MagnFieldHomY::MagnFieldHomY(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

 for (uint i=0;i<dim_in;i++) _refvalue[i]= qtymap_in.GetInputParser()->get("Bref");
}

//===========================================================================
MagnFieldHomZ::MagnFieldHomZ(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

 for (uint i=0;i<dim_in;i++) _refvalue[i]= qtymap_in.GetInputParser()->get("Bref");
}

//===========================================================================
MagnFieldHomAdjX::MagnFieldHomAdjX(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

  for (uint i=0;i<dim_in;i++) _refvalue[i]=1.;// qtymap_in.ml_prob.GetInputParser().get_par("Bref");
}

//===========================================================================
MagnFieldHomAdjY::MagnFieldHomAdjY(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

  for (uint i=0;i<dim_in;i++) _refvalue[i]=1.;// qtymap_in.ml_prob.GetInputParser().get_par("Bref");
}

//===========================================================================
MagnFieldHomAdjZ::MagnFieldHomAdjZ(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

  for (uint i=0;i<dim_in;i++) _refvalue[i]=1.;// qtymap_in.ml_prob.GetInputParser().get_par("Bref");
}

//==========================================================================
MagnFieldHomLagMult::MagnFieldHomLagMult(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 
  
  const double Bref = qtymap_in.GetInputParser()->get("Bref");
  const double Uref = qtymap_in.GetInputParser()->get("Uref");
  const double sigmaref = Uref*Bref;
  for (uint i=0;i<dim_in;i++) _refvalue[i]=sigmaref;
  
}

//==========================================================================
MagnFieldHomLagMultAdj::MagnFieldHomLagMultAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 
  
//   const double Bref = qtymap_in.ml_prob.GetInputParser().get_par("Bref");
//   const double Uref = qtymap_in.ml_prob.GetInputParser().get_par("Uref");
//   const double sigmaref = Uref*Bref;
  for (uint i=0;i<dim_in;i++) _refvalue[i]=1.;   //sigmaref;
  
}

//==========================================================================
MagnFieldExtX::MagnFieldExtX(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {

  for (uint i=0;i<dim_in;i++) _refvalue[i]= qtymap_in.GetInputParser()->get("Bref");
  
}

//==========================================================================
MagnFieldExtY::MagnFieldExtY(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {

  for (uint i=0;i<dim_in;i++) _refvalue[i]= qtymap_in.GetInputParser()->get("Bref");
  
}

//==========================================================================
MagnFieldExtZ::MagnFieldExtZ(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {

  for (uint i=0;i<dim_in;i++) _refvalue[i]= qtymap_in.GetInputParser()->get("Bref");
  
}

//===========================================================================
MagnFieldExtLagMult::MagnFieldExtLagMult(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {

  const double Bref = qtymap_in.GetInputParser()->get("Bref");
  const double Uref = qtymap_in.GetInputParser()->get("Uref");
  const double sigmaref = Uref*Bref;
  for (uint i=0;i<dim_in;i++) _refvalue[i]=sigmaref;

}

//===========================================================================
// an equation takes things from the equationsmap
//a quantity takes things from the quantity map

Pressure::Pressure(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

  
  for (uint i=0;i<dim_in;i++) _refvalue[i]= qtymap_in.GetInputParser()->get("pref");
}

//===========================================================================
PressureAdj::PressureAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

   for (uint i=0;i<dim_in;i++) _refvalue[i]= 1./*qtymap_in.ml_prob.GetInputParser()._pref*/;
}


//=========================================================================
VelocityX::VelocityX(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {  

    for (uint i=0;i<dim_in;i++) _refvalue[i] =  qtymap_in.GetInputParser()->get("Uref");
  
}

//=========================================================================
VelocityY::VelocityY(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {  

    for (uint i=0;i<dim_in;i++) _refvalue[i] =  qtymap_in.GetInputParser()->get("Uref");
  
}

//=========================================================================
VelocityZ::VelocityZ(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {  

    for (uint i=0;i<dim_in;i++) _refvalue[i] =  qtymap_in.GetInputParser()->get("Uref");
  
}

//=========================================================================
VelocityAdjX::VelocityAdjX(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {  

    for (uint i=0;i<dim_in;i++) _refvalue[i] = 1.;
    
}

//=========================================================================
VelocityAdjY::VelocityAdjY(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {  

    for (uint i=0;i<dim_in;i++) _refvalue[i] = 1.;
    
}

//=========================================================================
VelocityAdjZ::VelocityAdjZ(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {  

    for (uint i=0;i<dim_in;i++) _refvalue[i] = 1.;
    
}


//=================== END CONSTRUCTORS ================================
// ==================================================================
// ==================================================================



//===============================================
/// prescribed pressure at the boundary
//no initial condition for pressure is required, because it has no time derivative
//only the boundary condition in the Neumann part of the boundary has to be enforced
//so there is also no problem about consistency between IC and BC values,
// because there are NO IC values for pressure! 

void Pressure::Function_txyz(const double t, const double* xp,double* func) const {
  //the coordinates (x,y,z,t) of the VOLUME domain are NON-DIMENSIONAL
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

Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

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
  
 
  func[0] =  1./ _qtymap.GetInputParser()->get("pref")*( (le[1]-lb[1]) - /*x_rotshift*/xp[1] );

//this equation is in the reference frame CENTERED AT (0,0,0)  
  
  return;
  }
  
  
  

// ========================================================  
// ========================================================  
// ========================================================  


void VelocityX::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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


  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;    //u dot n
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll ) {  //right of the RefBox
    bc_flag[0]=0;    //u dot n
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )  {  //bottom  of the RefBox
     bc_flag[0]=0;      //u x n
   }
  
  if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  {  //top of the  of the RefBox
     bc_flag[0]=0;     //u x n
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
      bc_flag[0] = 0;   //u x n
  }
  
  if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  {
      bc_flag[0] = 0;      //u x n
  }
  
  
  
  
  return;
 
}



void VelocityY::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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


  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;    //u x n
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll ) {  //right of the RefBox
    bc_flag[0]=0;   //u x n
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )  {  //bottom  of the RefBox
      bc_flag[0]=0;   //u dot n   //leave this free for VELOCITY INLET
   }
  
  if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  {  //top of the  of the RefBox
//      bc_flag[0]=0;  //u dot n   //leave this free for outlet
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
      bc_flag[0] = 0;   //u x n
  }
  
  if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  {
      bc_flag[0] = 0;      //u x n
  }
  

  
  return;
 
}



void VelocityZ::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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


  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;    //u x n 
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll ) {  //right of the RefBox
    bc_flag[0]=0;   //u x n
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )  {  //bottom  of the RefBox
     bc_flag[0]=0;      //u x n
   }
  
  if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  {  //top of the  of the RefBox
     bc_flag[0]=0;     //u x n
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
      bc_flag[0]=0;                                               //u dot n  
  }
  
  if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  {
     bc_flag[0] = 0;                                                  //u dot n
  }
  
  return;
 
}



void Pressure::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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


  #if (DIMENSION==2)
  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
//   bc_flag[0]=0;
  }

 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
//  bc_flag[0]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
    bc_flag[0]=0;
                     //remember that this doesnt mean that the pressure at the boundary is fixed!
		     //its not a Dirichlet boundary condition for pressure!
		     //the initial guess of p at the boundary is changed by the equation after one step
		     //none of the pressure nodes are fixed, they are all COMPUTED after the first step
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
    bc_flag[0]=0;
  }
  

#elif (DIMENSION==3)

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
//      bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll ) {  //right of the RefBox
//      bc_flag[0]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )  {  //bottom  of the RefBox
//      bc_flag[0]=0;
  }
  
  if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  {  //top of the RefBox
      bc_flag[0]=0;     //tau dot n  //PRESSURE OUTLET
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
//      bc_flag[0]=0;
  }
  
  if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  {
//      bc_flag[0]=0;
  }

#endif
  
  return;
 
}



void MagnFieldHomX::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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


 #if (DIMENSION==2)

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {
//     bc_flag[0]=0;  //b.n useless with curl curl
  }
  
 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){
//    bc_flag[0]=0;  //b.n useless with curl curl
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  {
     bc_flag[0]=0;  //bxn
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {
    bc_flag[0]=0;  //bxn
 }
  
#else

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) { //  INSULATING
    bc_flag[0]=0;   //b.n
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){ //   INSULATING
    bc_flag[0]=0;   //b.n
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {
    bc_flag[0]=0;       //bxn
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {
    bc_flag[0]=0;         //bxn
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {   //  CONDUCTING, now insulating
   bc_flag[0] = 0;     //bxn
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {   //  CONDUCTING, now insulating
   bc_flag[0] = 0;     //bxn
  }
  
#endif
  
  
  return;
 
}


void MagnFieldHomY::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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


 #if (DIMENSION==2)

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {
     bc_flag[0]=0;  //bxn
  }
  
 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){
    bc_flag[0]=0;  //bxn
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  {
//      bc_flag[0]=0;      //leave this free
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {
//     bc_flag[0]=0;      //leave this free
 }
  
#else

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) { //  INSULATING
    bc_flag[0]=0;   //bxn
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){ //   INSULATING
    bc_flag[0]=0;   //bxn
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {
//     bc_flag[0]=0;    //b.n     //leave this free for inlet
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {
//     bc_flag[0]=0;      //b.n    //leave this free for outlet
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
   bc_flag[0] = 0;    //bxn 
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
   bc_flag[0] = 0;   //bxn
  }
  
#endif
  
  
  return;
 
}


void MagnFieldHomZ::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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


  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) { //  INSULATING
    bc_flag[0]=0;   //bxn
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){ //   INSULATING
    bc_flag[0]=0;   //bxn
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {
//     bc_flag[0]=0;    //bxn  //WHY ISNT THIS FIXED AS WELL?
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {
//     bc_flag[0]=0;     //bxn   //WHY ISNT THIS FIXED AS WELL?
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {   //  CONDUCTING, now insulating
//      bc_flag[0]=0;                                       //b.n 
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {   //  CONDUCTING, now insulating
//  bc_flag[0]=0;                           //b.n 
  }
  
  return;
 
}


void MagnFieldHomAdjX::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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
    bc_flag[0]=0;
  }
  
   if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
    bc_flag[0]=0;
  }
  
     if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox  
     bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox  
     bc_flag[0]=0;
  }
  
if ( x_rotshift[2] > -bdry_toll &&  x_rotshift[2] < bdry_toll ) { //current  
     bc_flag[0] = 0;
  }
  
  if ((le[2]-lb[2]) - x_rotshift[2] > -bdry_toll &&  (le[2]-lb[2]) -x_rotshift[2] < bdry_toll)  {
     bc_flag[0] = 0;
  }
  
  return;
 
}
 
 
void MagnFieldHomAdjY::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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
    bc_flag[0]=0;
  }
  
   if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
    bc_flag[0]=0;
  }
  
     if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox  
//      bc_flag[0]=0;    //u dot n 
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox  
//      bc_flag[0]=0;     //u dot n
  }
  
if ( x_rotshift[2] > -bdry_toll &&  x_rotshift[2] < bdry_toll ) { //current  
     bc_flag[0] = 0;
  }
  
  if ((le[2]-lb[2]) - x_rotshift[2] > -bdry_toll &&  (le[2]-lb[2]) -x_rotshift[2] < bdry_toll)  {
     bc_flag[0] = 0;
  }

  
  return;
 
}


void MagnFieldHomAdjZ::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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
    bc_flag[0]=0;
  }
  
   if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
    bc_flag[0]=0;
  }
  
     if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox  
     bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox  
     bc_flag[0]=0;
  }
  
if ( x_rotshift[2] > -bdry_toll &&  x_rotshift[2] < bdry_toll ) { //current  
//      bc_flag[0]=0;
  }
  
  if ((le[2]-lb[2]) - x_rotshift[2] > -bdry_toll &&  (le[2]-lb[2]) -x_rotshift[2] < bdry_toll)  {
//    bc_flag[0]=0;
  }
 
  
  return;
 
}

 
 
 void MagnFieldHomLagMult::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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


#if (DIMENSION==2)

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {
//   bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){
//   bc_flag[0]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  {
//   bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {
//   bc_flag[0]=0;
 }
  
#else

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {
//  bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){
//   bc_flag[0]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {
//   bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {
//   bc_flag[0]=0;
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
//     bc_flag[0]=0;
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
//     bc_flag[0]=0;
  }
  
#endif
  
  return;
 
}




void MagnFieldHomLagMultAdj::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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



  #if (DIMENSION==2)

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
//   bc_flag[0]=0;
  }
  
   if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
//   bc_flag[0]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox  
//   bc_flag[0]=0;      //NO BOUNDARY ADJOINT MHD PRESSURE
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox  
//   bc_flag[0]=0;      //NO BOUNDARY ADJOINT MHD PRESSURE 
  }
#else

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
//  bc_flag[0]=0;
  }
  
   if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
//  bc_flag[0]=0;
  }
  
     if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox  
//  bc_flag[0]=0;   //NO BOUNDARY ADJOINT MHD PRESSURE
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox  
//  bc_flag[0]=0;   //NO BOUNDARY ADJOINT MHD PRESSURE
  }
  
if ( x_rotshift[2] > -bdry_toll &&  x_rotshift[2] < bdry_toll ) { //current  
//      bc_flag[0]=0;
  }
  
  if ((le[2]-lb[2]) - x_rotshift[2] > -bdry_toll &&  (le[2]-lb[2]) -x_rotshift[2] < bdry_toll)  {
//      bc_flag[0]=0;
  }
#endif
  
  return;
 
}

 
 
 

void VelocityAdjX::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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


  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll) {  //right of the RefBox
    bc_flag[0]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll) {  //bottom  of the RefBox
     bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
     bc_flag[0]=0;
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
     bc_flag[0] = 0;   //u x n
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
    bc_flag[0] = 0;   //u x n
  }

  return;
 
}  
 
 
 

void VelocityAdjY::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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


  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll) {  //right of the RefBox
    bc_flag[0]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
//      bc_flag[0]=0;    //INSTEAD THIS MUST CORRESPOND TO THE DIRECT
                         //the SPACE of ADJOINT functions is the same as the SPACE for the DIRECT test functions
                         //if you fix this then you dont control...
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
//      bc_flag[0]=0;     //INSTEAD THIS MUST CORRESPOND TO THE DIRECT
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
     bc_flag[0] = 0;  //u x n
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
     bc_flag[0] = 0;  //u x n
  }


  return;
 
} 



 

void VelocityAdjZ::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
      bc_flag[0]=0;                                            //u dot n
  }
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
    bc_flag[0]=0;                                                //u dot n
  }


  return;
 
} 


 
 void PressureAdj::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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


#if (DIMENSION==2)
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
//   bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
//  bc_flag[0]=0;
  }

if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
//     bc_flag[0]=0;  //adjoint pressure: computing the integral is not correct, because the function p is prescribed
                    //at the boundary, so deltap = 0 at the boundary.
		    //You dont have to consider the symmetry with the direct equation!
                    //so THERE IS NO BOUNDARY INTEGRAL to be computed
		    //  COMMENT bc FOR THE ADJOINT PRESSURE
                    //The problem is that p_old is modified after one nonlinear step,
		    //so p was changing every time!
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
//    bc_flag[0]=0;   //  COMMENT bc FOR THE ADJOINT PRESSURE
  }
#else

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
//  bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox
//  bc_flag[0]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
//   bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
//   bc_flag[0]=0;
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
//   bc_flag[0]=0;
  }
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
//   bc_flag[0]=0;
  }
#endif
  
  
  return;
 
}  



void MagnFieldExtX::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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


  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox
    bc_flag[0]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
    bc_flag[0]=0;
   }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox
//     bc_flag[0] = 1;  //UseControl
  }
   if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
    bc_flag[0]=0;
  }
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
    bc_flag[0]=0;
  }
  
  return;
 
} 


//      bc_flag[1]=0;  //here when you control you must fix it,even if the corresponding b is let free, otherwise the control may use this as well
                      //well,wait,it depends: if you are using the Laplacian for MHD state, then Becont bc's must be consistent with the Laplacian
		      // given there.
		      //if you are using curlxcurl, then Becont should be consistent with curl-curl
		      //on the other hand, Be has a gamma*Laplacian FOR HERSELF... !
		      //so the Becont BC's should be consistent BOTH WITH MHD EQUATION and with MHDCONT EQUATION!
		      
		      
// wait, if you fix all dirichlet, then things go into Becontp pressure
//instead if you fix one Dir and one Neu then nothing goes into pressure

void MagnFieldExtY::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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



  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll)  {  //right of the RefBox
    bc_flag[0]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
    bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox
//     bc_flag[0] = 1; //UseControl;
  }
   if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
    bc_flag[0]=0;
  }
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
    bc_flag[0]=0;
  }
  
  
  return;
 
} 



void MagnFieldExtZ::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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


  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox
    bc_flag[0]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
    bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox
  
    bc_flag[0]=0;
  }
   if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
    bc_flag[0]=0;
  }
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
    bc_flag[0]=0;
  }
  
  
  return;
 
} 



void MagnFieldExtLagMult::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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

  
  
#if (DIMENSION==2)
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) { //left
//   bc_flag[0]=0;
  }
  
  if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){  //right
//  bc_flag[0]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
//   bc_flag[0]=0;     //comment it, don't do the integral
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {//top of the  of the RefBox
//   bc_flag[0]=0;       //comment it, don't do the integral

 }
  
#else

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
//     bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox
//     bc_flag[0]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
//     bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox
//     bc_flag[0]=0;
    
  }
   if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
//     bc_flag[0]=0;
  }
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
//     bc_flag[0]=0;
  }
  
#endif //DIMENSION
  
  
  return;
 
}

// =====================================================================
// ===================== INITIAL CONDITIONS ============================
// =====================================================================

void VelocityX::initialize_xyz(const double* xp, std::vector< double >& value) const {

    value[0] = 0.;
  
  return;
}





void VelocityY::initialize_xyz(const double* xp, std::vector< double >& value) const {

  const double Uref = _qtymap.GetInputParser()->get("Uref");
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
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

    value[0] = 0.;

   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  {  value[0] = (x_rotshift[0] - lb[0])*(le[0] - x_rotshift[0])*(x_rotshift[2] - lb[2])*(le[2]-x_rotshift[2])/Uref;  }
  
  
  return;
}



void VelocityZ::initialize_xyz(const double* xp, std::vector< double >& value) const {

    value[0] = 0.;
  
  return;
}





void Pressure::initialize_xyz(const double* xp, std::vector< double >& value) const {

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
  
#if (DIMENSION==2)
  
      Function_txyz(0.,&x_rotshift[0],&value[0]);
      
#elif (DIMENSION==3)
      
      value[0] = 0.;
  
#endif
  
  return;
}

void MagnFieldHomX::initialize_xyz(const double* xp, std::vector< double >& value) const {
  
  value[0] = 0.;
  return;
}


void MagnFieldHomY::initialize_xyz(const double* xp, std::vector< double >& value) const {
  
  value[0] = 0.;

  return;
}

void MagnFieldHomZ::initialize_xyz(const double* xp, std::vector< double >& value) const {
  
  value[0] = 0.;

  return;
}


void MagnFieldHomLagMult::initialize_xyz(const double* xp, std::vector< double >& value) const {
 
  value[0] = 0.;
  return;
}

void VelocityAdjX::initialize_xyz(const double* xp, std::vector< double >& value) const {

  value[0] = 0.;
  
  return;
}


void VelocityAdjY::initialize_xyz(const double* xp, std::vector< double >& value) const {

  value[0] = 0.;
  
  return;
}

void VelocityAdjZ::initialize_xyz(const double* xp, std::vector< double >& value) const {

  value[0] = 0.;
  
  return;
}

void PressureAdj::initialize_xyz(const double* xp, std::vector< double >& value) const {

  value[0] = 0.;
  
  return;
}

void MagnFieldHomAdjX::initialize_xyz(const double* xp, std::vector< double >& value) const {
  
  value[0] = 0.;
  return;
}


void MagnFieldHomAdjY::initialize_xyz(const double* xp, std::vector< double >& value) const {
  
  value[0] = 0.;
  return;
}

void MagnFieldHomAdjZ::initialize_xyz(const double* xp, std::vector< double >& value) const {
  
  value[0] = 0.;
  return;
}

void MagnFieldHomLagMultAdj::initialize_xyz(const double* xp, std::vector< double >& value) const {

  value[0] = 0.;

  return;
}

void MagnFieldExtX::initialize_xyz(const double* xp, std::vector< double >& value) const {

  const double Bref = _qtymap.GetInputParser()->get("Bref");
 
  value[0] = Bref/Bref;

  return;
}


void MagnFieldExtY::initialize_xyz(const double* xp, std::vector< double >& value) const {

  const double Bref = _qtymap.GetInputParser()->get("Bref");
 
  value[0] = 0./Bref;

  return;
}

void MagnFieldExtZ::initialize_xyz(const double* xp, std::vector< double >& value) const {

  const double Bref = _qtymap.GetInputParser()->get("Bref");
 
  value[0] = 0./Bref;

  return;
}


void MagnFieldExtLagMult::initialize_xyz(const double* xp, std::vector< double >& value) const {

  const double Uref = _qtymap.GetInputParser()->get("Uref");
  const double Bref = _qtymap.GetInputParser()->get("Bref");
 
  value[0] = 0./(Uref*Bref);

  return;
}
