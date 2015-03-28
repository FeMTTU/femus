
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
TempLift::TempLift(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in) 
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { }

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
VelocityX::VelocityX(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {  

   for (uint i=0;i<dim_in;i++) _refvalue[i] =  qtymap_in.GetInputParser()->get("Uref");
}

//=========================================================================
VelocityY::VelocityY(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {  

   for (uint i=0;i<dim_in;i++) _refvalue[i] =  qtymap_in.GetInputParser()->get("Uref");
}

//=================== END CONSTRUCTORS ================================
// ==================================================================
// ==================================================================


//=============================================================
void VelocityX::Function_txyz(const double /*t*/,const double* xp, double* func) const {

  
  func[0] = 0.;

  return;

}


//=============================================================
void VelocityY::Function_txyz(const double /*t*/,const double* xp, double* func) const {

  Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());
  
  const double rhof   = _qtymap.GetInputParser()->get("rho0");
  const double Uref   = _qtymap.GetInputParser()->get("Uref");
  const double Lref   = _qtymap.GetInputParser()->get("Lref");

  const double DpDz   = 1./*0.5*/;  //AAA: change it according to the pressure distribution!!!

  double DpDzad = DpDz*Lref/(rhof*Uref*Uref);

  const double magnitude = 5.*DpDzad*xp[0]*(1. - xp[0]);
 
  
  func[0] = magnitude;
  
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

Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());  //already nondimensionalized
 
  func[0] =  1./ _qtymap.GetInputParser()->get("pref")*( (box->_le[1] - box->_lb[1]) - xp[1] );

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

  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
     bc_flag[0]=0;
  }

 if ( (le[0]-lb[0]) - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll  ){ //right of the RefBox
    bc_flag[0]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
     bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox

//outflow only on part of the outlet
  if ( (x_rotshift[0]) > 0.70*(le[0]-lb[0]) ) {
      bc_flag[0]=0;      //ux
 }  //end part outflow
    else {  
      bc_flag[0]=0;  //ux
    }
  
  } //top RefBox
  

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

  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
     bc_flag[0]=0;
  }

 if ( (le[0]-lb[0]) - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll  ){ //right of the RefBox
    bc_flag[0]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
     bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox

//outflow only on part of the outlet
  if ( (x_rotshift[0]) > 0.70*(le[0]-lb[0]) ) {
//     bc_flag[0]=0;     //uy
 }  //end part outflow
    else {  
      bc_flag[0]=0;  //uy
    }
  
  } //top RefBox
  

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


  return;
 
}  
  
  
 void Temperature::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
// T' and its adjoint must be Dirichlet homogeneous everywhere on the boundary, by definition.


  const double bdry_toll = DEFAULT_BDRY_TOLL;
  

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


  const double bdry_toll = DEFAULT_BDRY_TOLL;
  

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
  
     bc_flag[0] = 0;
 
  }

 if ( (le[0]-lb[0])  - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll)  { //right of the RefBox
     bc_flag[0]=0;
   }

  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox

       if  ( (x_rotshift[0]) < 0.25*(le[0] - lb[0]) || ( x_rotshift[0]) > 0.75*(le[0] - lb[0]) )  {  bc_flag[0]=0; }
       
     }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox
      bc_flag[0]=0;
    }
 
  return;
}


 void TempAdj::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {

  const double bdry_toll = DEFAULT_BDRY_TOLL;
  

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
  
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;

  Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]); 

    value[0] = 1.;
  if ((box->_le[1]-box->_lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (box->_le[1]-box->_lb[1]) -(x_rotshift[1]) < bdry_toll)  {
    value[0] = 0.; 
  }


  return;
}

void VelocityX::initialize_xyz(const double* xp, std::vector< double >& value) const {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]); 
  
    value[0] = 0.;

    //==================================
    if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox

//below, inlet
  if ( (x_rotshift[0]) < 0.25*(box->_le[0] - box->_lb[0]) || ( x_rotshift[0]) > 0.75*(box->_le[0] - box->_lb[0]) ) { 
    value[0] = 0.;
  }
  else {
    value[0] = 0.; 
    }
  
  }
//============================================

//========================================
//left, inlet
 if  ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {
 
 if ( (x_rotshift[1]) > 0.4*(box->_le[1] - box->_lb[1]) && ( x_rotshift[1]) < 0.6*(box->_le[1]-box->_lb[1]) )  {  //left of the refbox
       value[0] = _qtymap.GetInputParser()->get("injsuc");    
      }
   }   
//============================================

//============================================
//====== outlet
  if ((box->_le[1]-box->_lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (box->_le[1]-box->_lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox

 if ( (x_rotshift[0]) < 0.71*(box->_le[0] - box->_lb[0]) ) {
   value[0] = 0.;
}

  }

  return;
}


void VelocityY::initialize_xyz(const double* xp, std::vector< double >& value) const {
  
  const double bdry_toll = DEFAULT_BDRY_TOLL;
  
  Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]); 
  
    value[0] = 0.;

 //==================================
    if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox

//below, inlet
  if ( (x_rotshift[0]) < 0.25*(box->_le[0] - box->_lb[0]) || ( x_rotshift[0]) > 0.75*(box->_le[0] - box->_lb[0]) ) { 
    value[0] = 0.;
  }
  else {
    value[0] = 1.; 
    }
  
  }
//============================================

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
