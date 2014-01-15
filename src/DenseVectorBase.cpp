// C++ includes
#include <iostream>
#include <iomanip> // for std::setw

// Local Includes
#include "DenseVectorBase.hpp"

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

