
//C++ includes
#include <cmath>
#include <iostream>

//library includes
#include "Typedefs.hpp"
#include "GeomEl.hpp"
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


// =================================================
void Temperature::Function_txyz(const double/* t*/, const double* xp,double* temp) const {

  const double Tref = _qtymap._physmap->get("Tref");

  Box* box = static_cast<Box*>(_qtymap._mesh.GetDomain());  
  
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

  return;
  }

 
  