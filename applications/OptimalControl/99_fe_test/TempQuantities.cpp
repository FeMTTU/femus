
//C++ includes
#include <cmath>
#include <iostream>

//library includes
#include "Typedefs.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "Box.hpp"

//application
#include "TempQuantities.hpp"

//=================== BEGIN CONSTRUCTORS ================================
// ==================================================================
// ==================================================================


//===========================================================================
Temperature::Temperature(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in) 
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {

}

//=================== END CONSTRUCTORS ================================
// ==================================================================
// ==================================================================



  
 
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
      bc_flag[0]=0;
  }

  if ( (le[0]-lb[0])  - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll)  { //right of the RefBox
      bc_flag[0]=0;
   }
   
  if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
//       bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox
//       bc_flag[0]=0;
  }

 if (_qtymap.GetMeshTwo()->get_dim() == 3) {
   
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
//      bc_flag[0]=0;
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
//       bc_flag[0]=0;
   }

}
  
  return;
}
 
 
  void Temperature::initialize_xyz(const double* xp, std::vector<double> & value) const {
   
    value[0] = 5.;
    
  }
  