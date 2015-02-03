#ifndef __physicsuser__
#define __physicsuser__

#include <iostream>
#include <cstdlib>

#include "Typedefs.hpp"
#include "Quantity.hpp"

using namespace femus;

class Temperature : public Quantity {

  public:
    
   Temperature(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~Temperature(){};
  
  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 
//specific function
  //this is the function of the IMPOSED DERIVATIVE of TEMPERATURE, aka heat flux
  void heatflux_txyz(const double t,const double* xyz, double* qflux) const;
  

};

#endif