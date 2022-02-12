#ifndef __femus_fe_test_TempQuantities_hpp__
#define __femus_fe_test_TempQuantities_hpp__

#include <iostream>
#include <cstdlib>

#include "Typedefs.hpp"
#include "Quantity.hpp"

using namespace femus;

class Temperature : public Quantity {

  public:
    
   Temperature(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};

#endif