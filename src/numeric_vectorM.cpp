#include "numeric_vectorM.hpp"

#include "FemusExtLib_conf.hpp"

#include "Typedefs_conf.hpp"

// Local Includes
#include "parallelM.hpp"
#include "petsc_vectorM.hpp"

// C++ includes
#include <cmath> // for std::abs
#include <memory>

//------------------------------------------------------------------
// NumericVectorM methods

// Full specialization for Real datatypes

std::auto_ptr<NumericVectorM >
NumericVectorM::build(const SolverPackage solver_package){
  // Build the appropriate vector
  switch (solver_package){
#ifdef FEMUS_HAVE_LASPACK
    case LASPACK_SOLVERSM:{
	std::auto_ptr<NumericVectorM > ap(new LaspackVectorM);
	return ap;
      }
#endif
#ifdef FEMUS_HAVE_PETSC
    case PETSC_SOLVERS:{
	std::auto_ptr<NumericVectorM > ap(new PetscVectorM);
	return ap;
      }
#endif
#ifdef LIBMESH_HAVE_TRILINOS
    case TRILINOS_SOLVERSM:{
	std::auto_ptr<NumericVectorM > ap(new EpetraVector<Real>);
	return ap;
      }
#endif
    default:
      std::auto_ptr<NumericVectorM > ap(NULL);
//       AutoPtr<NumericVectorM<T> > ap(new DistributedVector<T>);
      return ap;
    } 
  std::auto_ptr<NumericVectorM > ap(NULL);
  return ap;    
}



// Full specialization for float datatypes (DistributedVector wants this)

// template <>
// int NumericVectorM<float>::compare (const NumericVectorM<float> &other_vector,
// 				   const Real threshold) const
// {
//   assert (this->initialized());
//   assert (other_vector.initialized());
//   assert (this->first_local_index() == other_vector.first_local_index());
//   assert (this->last_local_index()  == other_vector.last_local_index());
// 
//   int rvalue     = -1;
//   unsigned int i = first_local_index();
// 
//   do
//     {
//       if ( std::abs( (*this)(i) - other_vector(i) ) > threshold )
// 	rvalue = i;
//       else
// 	i++;
//     }
//   while (rvalue==-1 && i<last_local_index());
// 
//   return rvalue;
// }

// Full specialization for double datatypes

int NumericVectorM::compare (const NumericVectorM &other_vector,
				    const Real threshold) const
{
  assert (this->initialized());
  assert (other_vector.initialized());
  assert (this->first_local_index() == other_vector.first_local_index());
  assert (this->last_local_index()  == other_vector.last_local_index());

  int rvalue     = -1;
  unsigned int i = first_local_index();

  do
    {
      if ( std::abs( (*this)(i) - other_vector(i) ) > threshold )
	rvalue = i;
      else
	i++;
    }
  while (rvalue==-1 && i<last_local_index());

  return rvalue;
}

#ifdef TRIPLE_PRECISION
// Full specialization for long double datatypes
template <>
int NumericVectorM<long double>::compare (const NumericVectorM<long double> &other_vector,
				         const Real threshold) const
{
  assert (this->initialized());
  assert (other_vector.initialized());
  assert (this->first_local_index() == other_vector.first_local_index());
  assert (this->last_local_index()  == other_vector.last_local_index());

  int rvalue     = -1;
  unsigned int i = first_local_index();

  do
    {
      if ( std::abs( (*this)(i) - other_vector(i) ) > threshold )
	rvalue = i;
      else
	i++;
    }
  while (rvalue==-1 && i<last_local_index());

  return rvalue;
}
#endif



// ---------------------------------------------------------
Real NumericVectorM::subset_l1_norm (const std::set<unsigned int> & indices){
  NumericVectorM & v = *this;
  
  std::set<unsigned int>::iterator it = indices.begin();
  const std::set<unsigned int>::iterator it_end = indices.end();
  Real norm = 0;
  for(; it!=it_end; ++it)    norm += std::abs(v(*it));
  ParallelM::sum(norm);
  return norm;
}

// --------------------------------------------------------
Real NumericVectorM::subset_l2_norm (const std::set<unsigned int> & indices){
  NumericVectorM & v = *this;
  std::set<unsigned int>::iterator it = indices.begin();
  const std::set<unsigned int>::iterator it_end = indices.end();
  Real norm = 0;
  for(; it!=it_end; ++it)    norm += (v(*it)*v(*it));
  ParallelM::sum(norm);
  return std::sqrt(norm);
}


Real NumericVectorM::subset_linfty_norm (const std::set<unsigned int> & indices){
  NumericVectorM & v = *this;
  std::set<unsigned int>::iterator it = indices.begin();
  const std::set<unsigned int>::iterator it_end = indices.end();
  Real norm = 0;
  for(; it!=it_end; ++it)    {
      Real value = std::abs(v(*it));
      if(value > norm)     norm = value;
    }
  ParallelM::max(norm); 

  return norm;
}



//------------------------------------------------------------------
// Explicit instantiations
// template class NumericVectorM<Number>;
