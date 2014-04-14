/*=========================================================================

 Program: FEMUS


} //end namespace femus


 Module: DenseVectorBase
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include "DenseVectorBase.hpp"


namespace femus {



// ===============================================================
void DenseVectorBase::print_scientific (std::ostream& os) const {
  // save the initial format flags
  std::ios_base::fmtflags os_flags = os.flags();
  // Print the vector entries.
  for (unsigned int i=0; i<this->size(); i++)
    os << std::setw(10)      << std::scientific       << std::setprecision(8)
       << this->el(i)       << std::endl;
  // reset the original format flags
  os.flags(os_flags);
}

// =================================================
void DenseVectorBase::print (std::ostream& os) const {
  for (unsigned int i=0; i<this->size(); i++)  os << std::setw(8) << this->el(i)  << std::endl;
}



} //end namespace femus


