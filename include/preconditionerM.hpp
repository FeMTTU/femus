#ifndef __preconditionerM_h__
#define __preconditionerM_h__

#include "FemusExtLib_conf.hpp"

#include "SolverPackage_enum.hpp"
#include "Precondtype_enum.hpp"


// C++ includes
#include <memory>
#include <iostream>
// Local includes
// #include "libmesh_common.h"
// #include "enum_solver_package.h"
// #include "enum_solver_type.h"
// #include "enum_preconditioner_type.h"
// #include "reference_counted_object.h"



// forward declarations
// template <typename T> class AutoPtr;
class SparseMatrixM;
class NumericVectorM;
// template <typename T> class ShellMatrix;





/**
 * This class provides a uniform interface for preconditioners.  This base
 * class is overloaded to provide linear solvers from different packages
 * like PETSC or Trilinos.
 *
 * In the below comments P is the matrix to be preconditioned with Apply()
 * performing the equivalent of the matrix vector product P^-1 x.  This
 * can also be thought of as (usually approximately) solving for Py=x.
 *
 * @author Derek Gaston, 2009
 */

class PreconditionerM 
 #ifdef LM_REFCOUNT
   : public ReferenceCountedObject<PreconditionerM>
 #endif
{
public:
  
  ///  Constructor. Initializes PreconditionerM data structures
  PreconditionerM ();
  /// Destructor.
  virtual ~PreconditionerM ();
  
  /// Builds a \p PreconditionerM using the linear solver package specified by \p solver_package
  static PreconditionerM * build(const SolverPackageM solver_package = PETSC_SOLVERSM);
  
  /// @returns true if the data structures 
  bool initialized () const { return _is_initialized; }

  /// Computes the preconditioned vector "y" based on input "x". Usually by solving Py=x to get the action of P^-1 x.
  virtual void apply(const NumericVectorM & x, NumericVectorM & y) = 0;
  
  /// Release all memory and clear data structures.
  virtual void clear () {}

  /// Initialize data structures if not done so already.
  virtual void init () {};

  /// Sets the matrix P to be preconditioned.
  void set_matrix(SparseMatrixM & mat);

  /// Returns the type of preconditioner to use.
  PreconditionerTypeM type () const{ return _preconditioner_type; }

  /// Sets the type of preconditioner to use.
  void set_type (const PreconditionerTypeM pct);
  
protected:

  /// The matrix P... ie the matrix to be preconditioned. This is often the actual system matrix of a linear sytem.
  SparseMatrixM * _matrix;
  
  /// Enum statitng with type of preconditioner to use.
  PreconditionerTypeM _preconditioner_type;
  
  /// Flag indicating if the data structures have been initialized.
  bool _is_initialized;
};




/*----------------------- inline functions ----------------------------------*/
// =============================================
inline PreconditionerM::PreconditionerM () :
  _matrix(NULL),
  _preconditioner_type (ILU_PRECONDM),
  _is_initialized      (false)
{
}
// =========================================================
inline PreconditionerM::~PreconditionerM (){  this->clear ();}
// ========================================================
inline void PreconditionerM::set_matrix(SparseMatrixM & mat){
  //If the matrix is changing then we (probably) need to reinitialize.
  _is_initialized = false;
  _matrix = &mat;
}
// ==========================================================
inline void PreconditionerM::set_type (const PreconditionerTypeM pct){
  //If the preconditioner type changes we (probably) need to reinitialize.
  _is_initialized = false;
  _preconditioner_type = pct;
}


#endif // #ifdef __preconditioner_h__
