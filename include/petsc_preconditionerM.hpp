#ifndef __petsc_preconditionerM_h__
#define __petsc_preconditionerM_h__


#include "FemusExtLib_conf.hpp"

#ifdef FEMUS_HAVE_PETSC

#include "PrecondtypeEnum.hpp"

// Local includes
#include "preconditionerM.hpp"

// Petsc includes
#include "petscpc.h"

// forward declarations
class SparseMatrixM;
class NumericVectorM;

/// This class provides an interface to  preconditioners  from Petsc.
class PetscPreconditionerM : public PreconditionerM{
  
  protected:
  PC _pc; ///<Preconditioner context
  Mat _mat;///< Petsc Matrix 
  
public:
  // Constructor
  ///  Constructor. Initializes PetscPreconditioner data structures 
  PetscPreconditionerM ():  PreconditionerM(){}
  /// Destructor.
  virtual ~PetscPreconditionerM (){this->clear ();  }  
  /// Release all memory and clear data structures. 
  virtual void clear () {}
  /// Initialize data structures if not done so already.
  virtual void init ();
  // Return
  /// Returns the actual Petsc PC struct.     
  PC pc() { return _pc; }
  
  // Compute 
   /// Computes the preconditioned vector "y" based on input "x".
  /// Usually by solving Py=x to get the action of P^-1 x.
  virtual void apply(const NumericVectorM & x, NumericVectorM & y);
  /// Tells PETSC to use the user-specified preconditioner  
  static void set_petsc_preconditioner_type 
             (const PreconditionerType & preconditioner_type, PC & pc);
};




/*----------------------- inline functions ----------------------------------*/
/*
inline PetscPreconditionerM::PetscPreconditionerM () :  PreconditionerM(){}
inline PetscPreconditionerM::~PetscPreconditionerM ()*/

#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // #ifdef __petsc_linear_solver_h__
