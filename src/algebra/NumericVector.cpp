/*=========================================================================

 Program: FEMUS
 Module: NumericVector
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
#include "NumericVector.hpp"
#include "FemusConfig.hpp"
#include "PetscVector.hpp"
#include "Parallel.hpp"
#include <cmath>
#include <memory>


namespace femus
{



//------------------------------------------------------------------
// NumericVector methods
  std::unique_ptr<NumericVector >
  NumericVector::build(const SolverPackage solver_package)
  {
    // Build the appropriate vector
    switch(solver_package) {
#ifdef HAVE_PETSC
      case PETSC_SOLVERS: {
          std::unique_ptr<NumericVector > ap(new PetscVector);
          return ap;
        }
#endif
#ifdef LIBMESH_HAVE_TRILINOS
      case TRILINOS_SOLVERS: {
          std::unique_ptr<NumericVector > ap(new EpetraVector<double>);
          return ap;
        }
#endif
      default:
        std::cerr << "SolverPackage solver_package:  Unrecognized: " << solver_package;
        abort();
    }
  }

//--------------------------------------------------------------------------------------------
  int NumericVector::compare(const NumericVector& other_vector,
                             const double threshold) const
  {
    assert(this->initialized());
    assert(other_vector.initialized());
    assert(this->first_local_index() == other_vector.first_local_index());
    assert(this->last_local_index()  == other_vector.last_local_index());

    int rvalue     = -1;
    int i = first_local_index();

    do {
      if(std::abs((*this)(i) - other_vector(i)) > threshold)
        rvalue = i;
      else
        i++;
    }
    while(rvalue == -1 && i < last_local_index());

    return rvalue;
  }

// ---------------------------------------------------------
  double NumericVector::subset_l1_norm(const std::set< int>& indices)
  {
    NumericVector& v = *this;

    std::set< int>::iterator it = indices.begin();
    const std::set< int>::iterator it_end = indices.end();
    double norm = 0;
    for(; it != it_end; ++it)    norm += std::abs(v(*it));
    Parallel::sum(norm);
    return norm;
  }

// --------------------------------------------------------
  double NumericVector::subset_l2_norm(const std::set< int>& indices)
  {
    NumericVector& v = *this;
    std::set< int>::iterator it = indices.begin();
    const std::set< int>::iterator it_end = indices.end();
    double norm = 0;
    for(; it != it_end; ++it)    norm += (v(*it) * v(*it));
    Parallel::sum(norm);
    return std::sqrt(norm);
  }

// --------------------------------------------------------
  double NumericVector::subset_linfty_norm(const std::set< int>& indices)
  {
    NumericVector& v = *this;
    std::set< int>::iterator it = indices.begin();
    const std::set< int>::iterator it_end = indices.end();
    double norm = 0;
    for(; it != it_end; ++it)    {
      double value = std::abs(v(*it));
      if(value > norm)     norm = value;
    }
    Parallel::max(norm);

    return norm;
  }

} //end namespace femus


