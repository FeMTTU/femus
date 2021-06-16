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
#include <limits>

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

    const int initialization_val = std::numeric_limits<int>::max();/*-1; this is not working in parallel*/

   std::cout << "---------------- " << initialization_val;
    int first_different_i     = initialization_val;
    int i = first_local_index();

    do {
      if( std::abs((*this)(i) - other_vector(i) ) > threshold)
        first_different_i = i;
      else
        i++;
    }
    while(first_different_i == initialization_val && i < last_local_index());

  // Find the correct first differing index in parallel
    Parallel::min(first_different_i);

//     if (n_procs > 1)
//     {
//       MPI_Allreduce (MPI_IN_PLACE, & first_different_i, 1,  MPI_INT, MPI_MIN,  MPI_COMM_WORLD);
//     }

  
  
  
  
  if (first_different_i == std::numeric_limits<int>::max())  { return -1; }

  return first_different_i;
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


